# =====================================================================
# BBT CONTINUOUS SMOOTH ROUTE v4.0
# 3D Slicer 5.12
#
# COMPLETE STANDALONE SCRIPT
#
# USER SELECTS
# ---------------------------------------------------------------------
# 1. Segmentation node
# 2. Airway segment
# 3. RouteStart markup
# 4. RouteTarget markup
# 5. Reference CT
#
# PROCESS
# ---------------------------------------------------------------------
# 1. Export selected airway
# 2. Identify connected components
# 3. Identify components containing RouteStart and RouteTarget
# 4. If disconnected, progressively bridge gaps
# 5. Store ALL bridges separately as Airway_INTERPOLATED
# 6. Create Airway_NAVIGATION_COMBINED
# 7. Verify Start -> Target connectivity
# 8. Skeletonise connected airway
# 9. Extract shortest Start -> Target skeleton path
# 10. Smooth path strongly to represent general airway direction
# 11. Constrain smoothed path to the Start-Target route component
# 12. Create Route_Centerline
#
# OUTPUTS
# ---------------------------------------------------------------------
# BBT_Navigation
#   ©À──Airway_INTERPOLATED
#   └──Airway_NAVIGATION_COMBINED
#
# Route_Centerline
#
# Gap_01_A
# Gap_01_B
# Gap_02_A
# Gap_02_B ...
#
# IMPORTANT
# ---------------------------------------------------------------------
# Airway_INTERPOLATED is NOT CT-derived anatomy.
#
# The source segmentation is NEVER modified.
#
# Route_Centerline is a SMOOTHED NAVIGATION TRAJECTORY.
# It is not represented as an exact anatomical centreline.
# =====================================================================


import slicer
import vtk
import qt

import numpy as np
import scipy.ndimage as ndi

import heapq
import math


# =====================================================================
# SETTINGS
# =====================================================================

OUTPUT_SEGMENTATION_NAME = "BBT_Navigation"

INTERPOLATED_SEGMENT_NAME = "Airway_INTERPOLATED"

COMBINED_SEGMENT_NAME = "Airway_NAVIGATION_COMBINED"

ROUTE_CURVE_NAME = "Route_Centerline"


# Maximum permitted individual artificial bridge
MAX_SINGLE_GAP_MM = 15.0

# Maximum number of bridges
MAX_BRIDGES = 10

# Bridge mask radius.
# Navigation connectivity only — NOT airway radius.
BRIDGE_RADIUS_MM = 0.75

BRIDGE_SAMPLES_PER_MM = 5.0


# =====================================================================
# MARKUP SNAPPING SETTINGS
# =====================================================================

# Maximum permitted distance from RouteStart / RouteTarget to the
# nearest voxel of the selected airway segment.
#
# Larger distances abort instead of silently snapping the route
# end to an unrelated airway voxel.
MAX_START_SNAP_MM = 5.0

MAX_TARGET_SNAP_MM = 5.0


# =====================================================================
# SMOOTHING SETTINGS
# =====================================================================

# Higher = smoother.
#
# Suggested:
# 5  = modest
# 10 = smooth
# 15 = very smooth
# 20 = aggressive
#
SMOOTHING_STRENGTH = 12.0

# Dense spline evaluation
DENSE_SAMPLE_SPACING_MM = 0.25

# Final Route_Centerline control point spacing
FINAL_POINT_SPACING_MM = 1.0


# =====================================================================
# HELPERS
# =====================================================================

def getNodeOrNone(name):

    try:
        return slicer.util.getNode(name)
    except:
        return None


def removeNodeIfExists(name):

    node = getNodeOrNone(name)

    if node:
        slicer.mrmlScene.RemoveNode(node)


def minimumDistanceMM(
    coordinatesKJI,
    pointKJI,
    spacingKJI
):
    """Smallest physical distance from voxel coordinates to a point."""

    delta = (
        np.asarray(coordinatesKJI, dtype=float)
        -
        np.asarray(pointKJI, dtype=float)
    ) * np.asarray(spacingKJI, dtype=float)

    return float(
        np.sqrt(
            np.min(
                np.sum(
                    delta * delta,
                    axis=1
                )
            )
        )
    )


def makesProgressTowardsTarget(
    candidateCoordinatesKJI,
    startTargetDistanceMM,
    targetKJI,
    spacingKJI
):
    """
    True when the candidate component reaches closer to the target
    than the start-connected component already does
    (startTargetDistanceMM).

    Components lying behind or beside the current route front are
    never bridged, so isolated fragments cannot consume bridges.
    """

    return (
        minimumDistanceMM(
            candidateCoordinatesKJI,
            targetKJI,
            spacingKJI
        )
        <
        startTargetDistanceMM
    )


def pointsOutsideMask(
    pointsKJI,
    mask
):
    """Number of points whose nearest voxel is outside the mask."""

    indices = np.round(
        np.asarray(pointsKJI, dtype=float)
    ).astype(int)

    outside = 0

    for k, j, i in indices:

        if not (
            0 <= k < mask.shape[0]
            and
            0 <= j < mask.shape[1]
            and
            0 <= i < mask.shape[2]
            and
            mask[k, j, i]
        ):
            outside += 1

    return outside


# =====================================================================
# ROBUST DROPDOWN
# =====================================================================

def chooseItem(title, text, items):

    if not items:
        raise RuntimeError(
            "No selectable items found."
        )

    dialog = qt.QDialog(
        slicer.util.mainWindow()
    )

    dialog.setWindowTitle(title)
    dialog.setMinimumWidth(520)

    layout = qt.QVBoxLayout(dialog)

    label = qt.QLabel(text)
    label.wordWrap = True
    layout.addWidget(label)

    combo = qt.QComboBox()

    for item in items:
        combo.addItem(str(item))

    layout.addWidget(combo)

    buttons = qt.QHBoxLayout()

    cancel = qt.QPushButton("Cancel")
    ok = qt.QPushButton("OK")

    buttons.addStretch(1)
    buttons.addWidget(cancel)
    buttons.addWidget(ok)

    layout.addLayout(buttons)

    state = {"accepted": False}

    def accept():

        state["accepted"] = True
        dialog.accept()

    def reject():

        dialog.reject()

    ok.connect(
        "clicked()",
        accept
    )

    cancel.connect(
        "clicked()",
        reject
    )

    dialog.exec()

    if not state["accepted"]:

        raise RuntimeError(
            "Selection cancelled."
        )

    return str(combo.currentText)


# =====================================================================
# YES / NO
# =====================================================================

def askYesNo(title, text, yesText="Continue"):

    box = qt.QMessageBox(
        slicer.util.mainWindow()
    )

    box.setWindowTitle(title)

    box.setText(text)

    box.setIcon(
        qt.QMessageBox.Warning
    )

    yes = box.addButton(
        yesText,
        qt.QMessageBox.AcceptRole
    )

    box.addButton(
        "Cancel",
        qt.QMessageBox.RejectRole
    )

    box.exec()

    return box.clickedButton() == yes


# =====================================================================
# COORDINATE CONVERSION
# =====================================================================

def rasToIJK(ras, volume):

    m = vtk.vtkMatrix4x4()

    volume.GetRASToIJKMatrix(m)

    ras4 = [
        float(ras[0]),
        float(ras[1]),
        float(ras[2]),
        1.0
    ]

    ijk4 = [
        0.0,
        0.0,
        0.0,
        1.0
    ]

    m.MultiplyPoint(
        ras4,
        ijk4
    )

    return np.array(
        ijk4[:3],
        dtype=float
    )


def ijkToRAS(ijk, volume):

    m = vtk.vtkMatrix4x4()

    volume.GetIJKToRASMatrix(m)

    ijk4 = [
        float(ijk[0]),
        float(ijk[1]),
        float(ijk[2]),
        1.0
    ]

    ras4 = [
        0.0,
        0.0,
        0.0,
        1.0
    ]

    m.MultiplyPoint(
        ijk4,
        ras4
    )

    return np.array(
        ras4[:3],
        dtype=float
    )


def rasToKJI(ras, volume):

    ijk = rasToIJK(
        ras,
        volume
    )

    return np.array(
        [
            ijk[2],
            ijk[1],
            ijk[0]
        ],
        dtype=float
    )


def kjiToRAS(kji, volume):

    return ijkToRAS(
        [
            kji[2],
            kji[1],
            kji[0]
        ],
        volume
    )


# =====================================================================
# SELECT SEGMENTATION
# =====================================================================

segNodes = slicer.util.getNodesByClass(
    "vtkMRMLSegmentationNode"
)

if not segNodes:

    raise RuntimeError(
        "No segmentation nodes found."
    )

choices = []
lookup = {}

for node in segNodes:

    text = (
        f"{node.GetName()} "
        f"[{node.GetID()}]"
    )

    choices.append(text)

    lookup[text] = node


selected = chooseItem(
    "BBT — Segmentation",
    "Select airway segmentation:",
    choices
)

segmentationNode = lookup[selected]


# =====================================================================
# SELECT AIRWAY SEGMENT
# =====================================================================

segmentation = (
    segmentationNode.GetSegmentation()
)

choices = []
lookup = {}

for i in range(
    segmentation.GetNumberOfSegments()
):

    sid = (
        segmentation.GetNthSegmentID(i)
    )

    segment = segmentation.GetSegment(
        sid
    )

    name = segment.GetName()

    text = (
        f"{name} [{sid}]"
    )

    choices.append(text)

    lookup[text] = (
        sid,
        name
    )


selected = chooseItem(
    "BBT — Airway",
    "Select airway segment:",
    choices
)

airwaySegmentID, airwaySegmentName = (
    lookup[selected]
)


# =====================================================================
# SELECT MARKUPS
# =====================================================================

markupNodes = []

for className in [
    "vtkMRMLMarkupsFiducialNode",
    "vtkMRMLMarkupsCurveNode"
]:

    nodes = slicer.util.getNodesByClass(
        className
    )

    for node in nodes:

        if node.GetNumberOfControlPoints() > 0:

            markupNodes.append(node)


if not markupNodes:

    raise RuntimeError(
        "No markup nodes found."
    )


choices = []
lookup = {}

for node in markupNodes:

    text = (
        f"{node.GetName()} "
        f"[{node.GetID()}]"
    )

    choices.append(text)

    lookup[text] = node


# =====================================================================
# START
# =====================================================================

selected = chooseItem(
    "BBT — RouteStart",
    "Select tracheal RouteStart:",
    choices
)

startNode = lookup[selected]


# =====================================================================
# TARGET
# =====================================================================

selected = chooseItem(
    "BBT — RouteTarget",
    "Select distal RouteTarget:",
    choices
)

targetNode = lookup[selected]


if startNode == targetNode:

    raise RuntimeError(
        "Start and Target cannot be the same node."
    )


startRAS = [0.0, 0.0, 0.0]
targetRAS = [0.0, 0.0, 0.0]

startNode.GetNthControlPointPositionWorld(
    0,
    startRAS
)

targetNode.GetNthControlPointPositionWorld(
    0,
    targetRAS
)

startRAS = np.array(
    startRAS,
    dtype=float
)

targetRAS = np.array(
    targetRAS,
    dtype=float
)


# =====================================================================
# SELECT CT
# =====================================================================

volumes = slicer.util.getNodesByClass(
    "vtkMRMLScalarVolumeNode"
)

if not volumes:

    raise RuntimeError(
        "No scalar volumes found."
    )

choices = []
lookup = {}

for node in volumes:

    text = (
        f"{node.GetName()} "
        f"[{node.GetID()}]"
    )

    choices.append(text)

    lookup[text] = node


selected = chooseItem(
    "BBT — Reference CT",
    "Select source CT:",
    choices
)

referenceVolume = lookup[selected]


if referenceVolume.GetParentTransformNode():

    raise RuntimeError(
        "Reference CT has a parent transform. "
        "Harden it first."
    )


# =====================================================================
# CONFIRM
# =====================================================================

text = (
    "CONTINUOUS SMOOTH BBT ROUTE\n\n"

    f"Segmentation:\n"
    f"{segmentationNode.GetName()}\n\n"

    f"Airway:\n"
    f"{airwaySegmentName}\n\n"

    f"Start:\n"
    f"{startNode.GetName()}\n\n"

    f"Target:\n"
    f"{targetNode.GetName()}\n\n"

    f"CT:\n"
    f"{referenceVolume.GetName()}\n\n"

    "The script will bridge disconnected airway fragments "
    "if required and then generate a smooth continuous "
    "Start-to-Target trajectory.\n\n"

    "All artificial connections remain explicitly separated "
    "from CT-derived anatomy."
)


if not askYesNo(
    "BBT",
    text,
    "Build route"
):

    raise RuntimeError(
        "Cancelled."
    )


# =====================================================================
# EXPORT AIRWAY
# =====================================================================

print("")
print("==============================================")
print("BBT CONTINUOUS SMOOTH ROUTE v4")
print("==============================================")


tempLabelmap = (
    slicer.mrmlScene.AddNewNodeByClass(
        "vtkMRMLLabelMapVolumeNode",
        "TEMP_BBT_Airway"
    )
)


ids = vtk.vtkStringArray()

ids.InsertNextValue(
    airwaySegmentID
)


success = (
    slicer.modules.segmentations.logic()
    .ExportSegmentsToLabelmapNode(
        segmentationNode,
        ids,
        tempLabelmap,
        referenceVolume
    )
)


if not success:

    raise RuntimeError(
        "Airway export failed."
    )


airway = (
    slicer.util.arrayFromVolume(
        tempLabelmap
    ) > 0
)


if not np.any(airway):

    raise RuntimeError(
        "Airway is empty."
    )


spacingIJK = np.array(
    referenceVolume.GetSpacing(),
    dtype=float
)

spacingKJI = spacingIJK[::-1]


# =====================================================================
# NEAREST AIRWAY VOXELS
# =====================================================================

airwayCoordinates = np.argwhere(
    airway
)


def nearestCoordinate(
    coordinates,
    ras
):

    target = rasToKJI(
        ras,
        referenceVolume
    )

    delta = (
        coordinates.astype(float)
        -
        target
    ) * spacingKJI

    d2 = np.sum(
        delta * delta,
        axis=1
    )

    index = int(
        np.argmin(d2)
    )

    return (
        tuple(
            coordinates[index]
        ),
        math.sqrt(
            float(d2[index])
        )
    )


startVoxel, startDistance = (
    nearestCoordinate(
        airwayCoordinates,
        startRAS
    )
)

targetVoxel, targetDistance = (
    nearestCoordinate(
        airwayCoordinates,
        targetRAS
    )
)


print("")
print(
    f"Start →airway: {startDistance:.2f} mm"
)

print(
    f"Target →airway: {targetDistance:.2f} mm"
)


if startDistance > MAX_START_SNAP_MM:

    raise RuntimeError(
        f"{startNode.GetName()} is {startDistance:.2f} mm "
        f"from the airway (limit {MAX_START_SNAP_MM:.2f} mm). "
        "Place it inside the tracheal lumen."
    )


if targetDistance > MAX_TARGET_SNAP_MM:

    raise RuntimeError(
        f"{targetNode.GetName()} is {targetDistance:.2f} mm "
        f"from the airway (limit {MAX_TARGET_SNAP_MM:.2f} mm). "
        "Place it inside the last CT-supported bronchus, "
        "not at the lesion."
    )


# =====================================================================
# WORKING MASK
#
# THIS VARIABLE IS NOW CREATED EXPLICITLY HERE.
# Therefore the previous NameError cannot occur.
# =====================================================================

workingMask = airway.copy()

bridgeMask = np.zeros_like(
    airway,
    dtype=bool
)

gapRecords = []

structure26 = np.ones(
    (3, 3, 3),
    dtype=np.uint8
)


# =====================================================================
# GAP BRIDGING LOOP
# =====================================================================

bridgeNumber = 0


while True:

    labels, numberComponents = ndi.label(
        workingMask,
        structure=structure26
    )

    startLabel = int(
        labels[startVoxel]
    )

    targetLabel = int(
        labels[targetVoxel]
    )

    print("")
    print(
        "Current components:",
        numberComponents
    )

    if startLabel == targetLabel:

        print(
            "Start and Target are connected."
        )

        break


    bridgeNumber += 1

    if bridgeNumber > MAX_BRIDGES:

        raise RuntimeError(
            "Maximum number of bridges reached."
        )


    startMask = (
        labels == startLabel
    )

    startTargetDistance = minimumDistanceMM(
        np.argwhere(
            startMask
        ),
        targetVoxel,
        spacingKJI
    )


    # ================================================================
    # DISTANCE FROM START-CONNECTED COMPONENT
    # ================================================================

    distanceToStart, nearestStartIndices = (
        ndi.distance_transform_edt(
            ~startMask,
            sampling=spacingKJI,
            return_indices=True
        )
    )


    # ================================================================
    # TARGET DIRECTION
    #
    # Used as a mild preference for components progressing generally
    # towards the selected target.
    # ================================================================

    startWorld = kjiToRAS(
        startVoxel,
        referenceVolume
    )

    targetWorld = kjiToRAS(
        targetVoxel,
        referenceVolume
    )

    routeVector = (
        targetWorld - startWorld
    )

    routeLength = np.linalg.norm(
        routeVector
    )

    if routeLength > 0:

        routeDirection = (
            routeVector / routeLength
        )

    else:

        routeDirection = np.array(
            [0.0, 0.0, 1.0]
        )


    bestScore = float("inf")

    bestDistance = None
    bestBKJI = None
    bestComponent = None


    componentIDs = np.unique(
        labels
    )


    for componentID in componentIDs:

        if (
            componentID == 0
            or
            componentID == startLabel
        ):
            continue


        componentMask = (
            labels == componentID
        )

        coords = np.argwhere(
            componentMask
        )

        if len(coords) == 0:
            continue


        if not makesProgressTowardsTarget(
            coords,
            startTargetDistance,
            targetVoxel,
            spacingKJI
        ):
            continue


        distances = (
            distanceToStart[
                componentMask
            ]
        )

        localIndex = int(
            np.argmin(distances)
        )

        gapDistance = float(
            distances[localIndex]
        )

        candidateKJI = coords[
            localIndex
        ]

        candidateRAS = kjiToRAS(
            candidateKJI,
            referenceVolume
        )


        # ============================================================
        # TARGET PROGRESSION PENALTY
        #
        # This is deliberately mild.
        # Distance remains the dominant criterion.
        # ============================================================

        candidateVector = (
            candidateRAS - startWorld
        )

        candidateLength = np.linalg.norm(
            candidateVector
        )


        directionPenalty = 0.0


        if candidateLength > 0:

            candidateDirection = (
                candidateVector
                /
                candidateLength
            )

            dot = float(
                np.clip(
                    np.dot(
                        candidateDirection,
                        routeDirection
                    ),
                    -1.0,
                    1.0
                )
            )

            # 0 for same direction;
            # increasing penalty for backwards direction.

            directionPenalty = (
                (1.0 - dot) * 2.0
            )


        score = (
            gapDistance
            +
            directionPenalty
        )


        if score < bestScore:

            bestScore = score
            bestDistance = gapDistance
            bestBKJI = candidateKJI
            bestComponent = int(
                componentID
            )


    if bestComponent is None:

        raise RuntimeError(
            "No disconnected component lies closer to the target "
            "than the start-connected airway."
        )


    gapBKJI = np.array(
        bestBKJI,
        dtype=int
    )


    gapAKJI = np.array(
        [
            nearestStartIndices[0][tuple(gapBKJI)],
            nearestStartIndices[1][tuple(gapBKJI)],
            nearestStartIndices[2][tuple(gapBKJI)]
        ],
        dtype=int
    )


    gapARAS = kjiToRAS(
        gapAKJI,
        referenceVolume
    )

    gapBRAS = kjiToRAS(
        gapBKJI,
        referenceVolume
    )


    gapLength = float(
        bestDistance
    )


    print(
        f"Candidate gap {bridgeNumber}: "
        f"{gapLength:.2f} mm"
    )


    # ================================================================
    # GAP LIMIT
    # ================================================================

    if gapLength > MAX_SINGLE_GAP_MM:

        raise RuntimeError(
            f"Required gap is {gapLength:.2f} mm, "
            f"above limit {MAX_SINGLE_GAP_MM:.2f} mm."
        )


    # ================================================================
    # AUTHORISATION
    #
    # Centre the slice views on the gap so the CT can be inspected
    # before deciding.
    # ================================================================

    gapMidRAS = (
        gapARAS + gapBRAS
    ) / 2.0

    try:
        slicer.modules.markups.logic().JumpSlicesToLocation(
            float(gapMidRAS[0]),
            float(gapMidRAS[1]),
            float(gapMidRAS[2]),
            True
        )
        slicer.app.processEvents()
    except Exception as error:
        print(
            "WARNING: could not centre slice views on gap:",
            error
        )

    text = (
        f"GAP {bridgeNumber:02d}\n\n"

        f"Gap length: {gapLength:.2f} mm\n\n"

        "Slice views are centred on the gap midpoint "
        f"(RAS {gapMidRAS[0]:.1f}, {gapMidRAS[1]:.1f}, "
        f"{gapMidRAS[2]:.1f}). Check the CT before deciding.\n\n"

        "This geometric bridge is NOT CT-derived anatomy.\n\n"

        "It will be stored separately in "
        "Airway_INTERPOLATED.\n\n"

        "Authorise this bridge?"
    )


    if not askYesNo(
        "BBT — Interpolation",
        text,
        "Authorise"
    ):

        raise RuntimeError(
            "Interpolation declined."
        )


    # ================================================================
    # GAP MARKERS
    # ================================================================

    gapAName = (
        f"Gap_{bridgeNumber:02d}_A"
    )

    gapBName = (
        f"Gap_{bridgeNumber:02d}_B"
    )


    removeNodeIfExists(
        gapAName
    )

    removeNodeIfExists(
        gapBName
    )


    gapA = (
        slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLMarkupsFiducialNode",
            gapAName
        )
    )

    gapA.AddControlPointWorld(
        vtk.vtkVector3d(
            *[
                float(x)
                for x in gapARAS
            ]
        )
    )


    gapB = (
        slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLMarkupsFiducialNode",
            gapBName
        )
    )

    gapB.AddControlPointWorld(
        vtk.vtkVector3d(
            *[
                float(x)
                for x in gapBRAS
            ]
        )
    )


    # ================================================================
    # BUILD BRIDGE CENTRE
    # ================================================================

    bridgeCentre = np.zeros_like(
        airway,
        dtype=bool
    )


    samples = max(
        2,
        int(
            math.ceil(
                gapLength
                *
                BRIDGE_SAMPLES_PER_MM
            )
        )
    )


    for fraction in np.linspace(
        0.0,
        1.0,
        samples
    ):

        ras = (
            gapARAS * (1.0 - fraction)
            +
            gapBRAS * fraction
        )

        kji = np.round(
            rasToKJI(
                ras,
                referenceVolume
            )
        ).astype(int)

        k, j, i = kji

        if (
            0 <= k < airway.shape[0]
            and
            0 <= j < airway.shape[1]
            and
            0 <= i < airway.shape[2]
        ):

            bridgeCentre[
                k,
                j,
                i
            ] = True


    # ================================================================
    # BRIDGE THICKNESS
    # ================================================================

    bridgeDistance = (
        ndi.distance_transform_edt(
            ~bridgeCentre,
            sampling=spacingKJI
        )
    )


    thisBridge = (
        bridgeDistance
        <=
        BRIDGE_RADIUS_MM
    )


    thisBridge[
        tuple(gapAKJI)
    ] = True

    thisBridge[
        tuple(gapBKJI)
    ] = True


    bridgeMask = np.logical_or(
        bridgeMask,
        thisBridge
    )


    # ================================================================
    # UPDATE WORKING MASK
    # ================================================================

    workingMask = np.logical_or(
        airway,
        bridgeMask
    )


    gapRecords.append(
        {
            "number": bridgeNumber,
            "length": gapLength,
            "A": gapARAS.copy(),
            "B": gapBRAS.copy()
        }
    )


# =====================================================================
# VERIFY FINAL CONNECTIVITY
# =====================================================================

labels, _ = ndi.label(
    workingMask,
    structure=structure26
)


if (
    labels[startVoxel]
    !=
    labels[targetVoxel]
):

    raise RuntimeError(
        "Start and Target remain disconnected."
    )


# =====================================================================
# ISOLATE ROUTE COMPONENT
# =====================================================================

routeLabel = int(
    labels[startVoxel]
)


routeComponent = (
    labels == routeLabel
)


# =====================================================================
# SKELETONISE
# =====================================================================

print("")
print(
    "Skeletonising continuous airway..."
)


try:

    from skimage.morphology import skeletonize

except ImportError:

    slicer.util.pip_install(
        "scikit-image"
    )

    from skimage.morphology import skeletonize


skeleton = np.asarray(
    skeletonize(
        routeComponent
    ),
    dtype=bool
)


skeletonCoordinates = np.argwhere(
    skeleton
)


if len(skeletonCoordinates) == 0:

    raise RuntimeError(
        "Skeletonisation failed."
    )


# =====================================================================
# NEAREST SKELETON VOXELS
# =====================================================================

def nearestSkeleton(ras):

    target = rasToKJI(
        ras,
        referenceVolume
    )

    delta = (
        skeletonCoordinates.astype(float)
        -
        target
    ) * spacingKJI

    d2 = np.sum(
        delta * delta,
        axis=1
    )

    return tuple(
        skeletonCoordinates[
            int(np.argmin(d2))
        ]
    )


skeletonStart = nearestSkeleton(
    startRAS
)

skeletonTarget = nearestSkeleton(
    targetRAS
)


# =====================================================================
# BUILD GRAPH
# =====================================================================

neighbours = []


for dk in [-1, 0, 1]:

    for dj in [-1, 0, 1]:

        for di in [-1, 0, 1]:

            if (
                dk == 0
                and
                dj == 0
                and
                di == 0
            ):
                continue

            cost = math.sqrt(
                (dk * spacingKJI[0]) ** 2
                +
                (dj * spacingKJI[1]) ** 2
                +
                (di * spacingKJI[2]) ** 2
            )

            neighbours.append(
                (
                    dk,
                    dj,
                    di,
                    cost
                )
            )


def heuristic(a, b):

    delta = (
        np.array(a, dtype=float)
        -
        np.array(b, dtype=float)
    ) * spacingKJI

    return float(
        np.linalg.norm(delta)
    )


# =====================================================================
# A* ROUTE
# =====================================================================

print(
    "Calculating Start →Target route..."
)


openHeap = [
    (
        heuristic(
            skeletonStart,
            skeletonTarget
        ),
        0.0,
        skeletonStart
    )
]


cameFrom = {}

gScore = {
    skeletonStart: 0.0
}

visited = set()

shape = skeleton.shape

routeFound = False


while openHeap:

    _, currentG, current = (
        heapq.heappop(
            openHeap
        )
    )

    if current in visited:
        continue

    visited.add(current)

    if current == skeletonTarget:

        routeFound = True
        break


    ck, cj, ci = current


    for (
        dk,
        dj,
        di,
        stepCost
    ) in neighbours:

        nk = ck + dk
        nj = cj + dj
        ni = ci + di

        if (
            nk < 0
            or nj < 0
            or ni < 0
            or nk >= shape[0]
            or nj >= shape[1]
            or ni >= shape[2]
        ):
            continue

        if not skeleton[
            nk,
            nj,
            ni
        ]:
            continue

        neighbour = (
            nk,
            nj,
            ni
        )

        tentative = (
            currentG
            +
            stepCost
        )

        if tentative < gScore.get(
            neighbour,
            float("inf")
        ):

            cameFrom[
                neighbour
            ] = current

            gScore[
                neighbour
            ] = tentative

            score = (
                tentative
                +
                heuristic(
                    neighbour,
                    skeletonTarget
                )
            )

            heapq.heappush(
                openHeap,
                (
                    score,
                    tentative,
                    neighbour
                )
            )


if not routeFound:

    raise RuntimeError(
        "Could not extract continuous skeleton route."
    )


# =====================================================================
# RECONSTRUCT RAW ROUTE
# =====================================================================

routeVoxels = [
    skeletonTarget
]

current = skeletonTarget


while current != skeletonStart:

    current = cameFrom[
        current
    ]

    routeVoxels.append(
        current
    )


routeVoxels.reverse()


routeRAS = np.array(
    [
        kjiToRAS(
            voxel,
            referenceVolume
        )
        for voxel in routeVoxels
    ],
    dtype=float
)


# =====================================================================
# CLEAN DUPLICATE POINTS
# =====================================================================

cleanRoute = [
    routeRAS[0]
]


for point in routeRAS[1:]:

    if (
        np.linalg.norm(
            point - cleanRoute[-1]
        )
        >
        0.05
    ):

        cleanRoute.append(
            point
        )


cleanRoute = np.asarray(
    cleanRoute,
    dtype=float
)


if len(cleanRoute) < 4:

    raise RuntimeError(
        "Insufficient points for smoothing."
    )


# =====================================================================
# SMOOTH SPLINE
# =====================================================================

print("")
print(
    "Creating smooth navigation trajectory..."
)


from scipy.interpolate import splprep, splev


segmentLengths = np.linalg.norm(
    np.diff(
        cleanRoute,
        axis=0
    ),
    axis=1
)


cumulative = np.concatenate(
    (
        [0.0],
        np.cumsum(
            segmentLengths
        )
    )
)


rawLength = float(
    cumulative[-1]
)


u = (
    cumulative
    /
    rawLength
)


smoothingParameter = (
    SMOOTHING_STRENGTH
    *
    len(cleanRoute)
    *
    0.1
)


tck, _ = splprep(
    [
        cleanRoute[:, 0],
        cleanRoute[:, 1],
        cleanRoute[:, 2]
    ],
    u=u,
    s=smoothingParameter,
    k=3
)


denseCount = max(
    20,
    int(
        math.ceil(
            rawLength
            /
            DENSE_SAMPLE_SPACING_MM
        )
    )
)


denseU = np.linspace(
    0.0,
    1.0,
    denseCount
)


xyz = splev(
    denseU,
    tck
)


smoothRoute = np.column_stack(
    xyz
)


smoothRoute[0] = cleanRoute[0]
smoothRoute[-1] = cleanRoute[-1]


# =====================================================================
# CONSTRAIN SPLINE TO ROUTE COMPONENT
#
# Only the Start-to-Target connected component is used. Constraining
# to the whole navigation mask could project points onto a
# neighbouring branch or an isolated fragment.
#
# Projection moves a point to the nearest interior voxel, which lies
# at the airway wall rather than the lumen centre.
# =====================================================================

print(
    "Constraining smooth route to route component..."
)


outsideDistance, nearestInside = (
    ndi.distance_transform_edt(
        ~routeComponent,
        sampling=spacingKJI,
        return_indices=True
    )
)


constrainedRoute = []

projectedPoints = 0


for point in smoothRoute:

    kji = np.round(
        rasToKJI(
            point,
            referenceVolume
        )
    ).astype(int)


    k = int(
        np.clip(
            kji[0],
            0,
            routeComponent.shape[0] - 1
        )
    )

    j = int(
        np.clip(
            kji[1],
            0,
            routeComponent.shape[1] - 1
        )
    )

    i = int(
        np.clip(
            kji[2],
            0,
            routeComponent.shape[2] - 1
        )
    )


    if routeComponent[
        k,
        j,
        i
    ]:

        constrainedRoute.append(
            point
        )

    else:

        nk = int(
            nearestInside[0][
                k,
                j,
                i
            ]
        )

        nj = int(
            nearestInside[1][
                k,
                j,
                i
            ]
        )

        ni = int(
            nearestInside[2][
                k,
                j,
                i
            ]
        )


        constrainedRoute.append(
            kjiToRAS(
                (
                    nk,
                    nj,
                    ni
                ),
                referenceVolume
            )
        )

        projectedPoints += 1


constrainedRoute = np.asarray(
    constrainedRoute,
    dtype=float
)


# =====================================================================
# MOVING-AVERAGE SMOOTHING AFTER CONSTRAINT
# =====================================================================

if len(constrainedRoute) >= 7:

    filtered = constrainedRoute.copy()

    weights = np.array(
        [
            1,
            2,
            3,
            4,
            3,
            2,
            1
        ],
        dtype=float
    )

    weights /= np.sum(weights)


    for index in range(
        3,
        len(constrainedRoute) - 3
    ):

        neighbourhood = (
            constrainedRoute[
                index - 3:
                index + 4
            ]
        )

        candidate = np.sum(
            neighbourhood
            *
            weights[:, None],
            axis=0
        )


        kji = np.round(
            rasToKJI(
                candidate,
                referenceVolume
            )
        ).astype(int)


        k, j, i = kji


        if (
            0 <= k < routeComponent.shape[0]
            and
            0 <= j < routeComponent.shape[1]
            and
            0 <= i < routeComponent.shape[2]
            and
            routeComponent[
                k,
                j,
                i
            ]
        ):

            filtered[index] = candidate


    constrainedRoute = filtered


# =====================================================================
# RESAMPLE BY DISTANCE
# =====================================================================

lengths = np.linalg.norm(
    np.diff(
        constrainedRoute,
        axis=0
    ),
    axis=1
)


cum = np.concatenate(
    (
        [0.0],
        np.cumsum(lengths)
    )
)


finalLength = float(
    cum[-1]
)


sampleDistances = np.arange(
    0.0,
    finalLength,
    FINAL_POINT_SPACING_MM
)


sampleDistances = np.append(
    sampleDistances,
    finalLength
)


finalRoute = []


for distance in sampleDistances:

    index = np.searchsorted(
        cum,
        distance
    )


    if index == 0:

        finalRoute.append(
            constrainedRoute[0]
        )

        continue


    if index >= len(
        constrainedRoute
    ):

        finalRoute.append(
            constrainedRoute[-1]
        )

        continue


    d0 = cum[index - 1]
    d1 = cum[index]


    if d1 == d0:

        fraction = 0.0

    else:

        fraction = (
            (distance - d0)
            /
            (d1 - d0)
        )


    point = (
        constrainedRoute[index - 1]
        *
        (1.0 - fraction)
        +
        constrainedRoute[index]
        *
        fraction
    )


    finalRoute.append(
        point
    )


finalRoute = np.asarray(
    finalRoute,
    dtype=float
)


# =====================================================================
# VERIFY FINAL ROUTE AGAINST ROUTE COMPONENT
#
# Resampled points are linear interpolations between constrained
# points and are not guaranteed to remain inside the airway.
# =====================================================================

finalPointsOutside = pointsOutsideMask(
    [
        rasToKJI(
            point,
            referenceVolume
        )
        for point in finalRoute
    ],
    routeComponent
)

if finalPointsOutside:

    print(
        "WARNING:",
        finalPointsOutside,
        "of",
        len(finalRoute),
        "final route points lie outside the route component. "
        "Review Route_Centerline against the CT."
    )


# =====================================================================
# CREATE ROUTE CENTERLINE
# =====================================================================

removeNodeIfExists(
    ROUTE_CURVE_NAME
)


routeCurve = (
    slicer.mrmlScene.AddNewNodeByClass(
        "vtkMRMLMarkupsCurveNode",
        ROUTE_CURVE_NAME
    )
)


for point in finalRoute:

    routeCurve.AddControlPointWorld(
        vtk.vtkVector3d(
            float(point[0]),
            float(point[1]),
            float(point[2])
        )
    )


# Coordinates are already smoothed.
# Keep Slicer interpolation linear so another spline does not cut
# through airway walls.

try:

    routeCurve.SetCurveTypeToLinear()

except:

    pass


# =====================================================================
# CREATE OUTPUT SEGMENTATION
# =====================================================================

removeNodeIfExists(
    OUTPUT_SEGMENTATION_NAME
)


outputSegmentation = (
    slicer.mrmlScene.AddNewNodeByClass(
        "vtkMRMLSegmentationNode",
        OUTPUT_SEGMENTATION_NAME
    )
)


outputSegmentation.CreateDefaultDisplayNodes()


# =====================================================================
# ARRAY -> SEGMENT
# =====================================================================

def importArrayAsSegment(
    array,
    segmentName
):

    temp = (
        slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLLabelMapVolumeNode",
            "TEMP_" + segmentName
        )
    )


    matrix = vtk.vtkMatrix4x4()

    referenceVolume.GetIJKToRASMatrix(
        matrix
    )

    temp.SetIJKToRASMatrix(
        matrix
    )


    slicer.util.updateVolumeFromArray(
        temp,
        array.astype(
            np.uint8
        )
    )


    outputSegments = (
        outputSegmentation.GetSegmentation()
    )


    before = set(
        outputSegments.GetNthSegmentID(i)
        for i in range(
            outputSegments.GetNumberOfSegments()
        )
    )


    slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(
        temp,
        outputSegmentation
    )


    newIDs = [
        outputSegments.GetNthSegmentID(i)
        for i in range(
            outputSegments.GetNumberOfSegments()
        )
        if outputSegments.GetNthSegmentID(i)
        not in before
    ]


    slicer.mrmlScene.RemoveNode(
        temp
    )


    if len(newIDs) != 1:

        raise RuntimeError(
            "Failed creating "
            +
            segmentName
        )


    sid = newIDs[0]


    outputSegments.GetSegment(
        sid
    ).SetName(
        segmentName
    )


    return sid


# =====================================================================
# INTERPOLATION OUTPUT
# =====================================================================

if np.any(
    bridgeMask
):

    importArrayAsSegment(
        bridgeMask,
        INTERPOLATED_SEGMENT_NAME
    )


# =====================================================================
# COMBINED OUTPUT
# =====================================================================

importArrayAsSegment(
    workingMask,
    COMBINED_SEGMENT_NAME
)


outputSegmentation.CreateClosedSurfaceRepresentation()


# =====================================================================
# AUDIT INFORMATION
# =====================================================================

routeCurve.SetAttribute(
    "BBT.Type",
    "Smoothed navigation trajectory"
)

routeCurve.SetAttribute(
    "BBT.SourceSegmentation",
    segmentationNode.GetName()
)

routeCurve.SetAttribute(
    "BBT.SourceAirwaySegment",
    airwaySegmentName
)

routeCurve.SetAttribute(
    "BBT.ReferenceCT",
    referenceVolume.GetName()
)

routeCurve.SetAttribute(
    "BBT.InterpolationUsed",
    str(
        len(gapRecords) > 0
    )
)

routeCurve.SetAttribute(
    "BBT.NumberOfInterpolatedGaps",
    str(
        len(gapRecords)
    )
)

routeCurve.SetAttribute(
    "BBT.SmoothingStrength",
    str(
        SMOOTHING_STRENGTH
    )
)

routeCurve.SetAttribute(
    "BBT.Constraint",
    "Airway_NAVIGATION_COMBINED Start-Target component"
)

routeCurve.SetAttribute(
    "BBT.PointsOutsideRouteComponent",
    str(
        finalPointsOutside
    )
)

# Control points run from RouteStart to RouteTarget.
# The virtual bronchoscopy script uses this instead of guessing.
routeCurve.SetAttribute(
    "BBT.RouteDirection",
    "StartToTarget"
)


outputSegmentation.SetAttribute(
    "BBT.InterpolationIsCTDerived",
    "False"
)

outputSegmentation.SetAttribute(
    "BBT.NumberOfBridges",
    str(
        len(gapRecords)
    )
)


# =====================================================================
# CLEAN TEMPORARY NODE
# =====================================================================

slicer.mrmlScene.RemoveNode(
    tempLabelmap
)


# =====================================================================
# FINAL REPORT
# =====================================================================

print("")
print("==============================================")
print("BBT ROUTE COMPLETE")
print("==============================================")
print("")

print(
    "Interpolated gaps:",
    len(gapRecords)
)


for record in gapRecords:

    print(
        f"Gap {record['number']:02d}: "
        f"{record['length']:.2f} mm"
    )


print(
    f"Raw skeleton route: "
    f"{rawLength:.1f} mm"
)

print(
    f"Smooth trajectory: "
    f"{finalLength:.1f} mm"
)

print(
    "Final route points:",
    len(finalRoute)
)

print(
    "Spline points projected back into airway:",
    projectedPoints
)

print(
    "Final points outside route component:",
    finalPointsOutside
)


# =====================================================================
# FINAL MESSAGE
# =====================================================================

gapText = ""


if gapRecords:

    gapText = "\n".join(
        [
            f"Gap {r['number']:02d}: "
            f"{r['length']:.2f} mm"
            for r in gapRecords
        ]
    )

else:

    gapText = "No interpolation required."


qt.QMessageBox.information(
    slicer.util.mainWindow(),
    "BBT — Complete",
    (
        "Continuous smooth route created.\n\n"

        f"{gapText}\n\n"

        f"Route length: "
        f"{finalLength:.1f} mm\n\n"

        "Created:\n"
        "• Route_Centerline\n"
        "• Airway_NAVIGATION_COMBINED\n"
        +
        (
            "• Airway_INTERPOLATED\n"
            if gapRecords
            else ""
        )
        +
        "\nRoute_Centerline is a smoothed navigation "
        "trajectory representing the general airway direction.\n\n"

        "Interpolated regions are NOT CT-derived anatomy."
    )
)
