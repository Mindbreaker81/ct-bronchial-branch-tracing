# =====================================================================
# BBT VIRTUAL BRONCHOSCOPY v7.0
# PATIENT-FIXED ANATOMICAL CAMERA
# 3D Slicer 5.12
#
# REQUIRED:
#   Route_Centerline
#
# OPTIONAL:
#   B1, B2, B3... Markups Fiducials
#   Description = desired clock position (1-12)
#
# CAMERA MODEL
# ---------------------------------------------------------------------
# Patient remains FIXED in CT/RAS anatomical coordinates.
#
# Slicer RAS:
#   +X = Right
#   +Y = Anterior
#   +Z = Superior / Head
#
# Therefore:
#   Head -> Feet = -Z
#
# The patient/model NEVER rotates.
#
# Camera:
#   1. Starts at tracheal end of Route_Centerline.
#   2. Route direction is taken from BBT.RouteDirection when script 03
#      set it; otherwise it is chosen so that initial travel is
#      preferentially HEAD -> FEET.
#   3. Forward vector follows Route_Centerline.
#   4. Anatomical anterior defines initial camera "up".
#   5. Camera frame is parallel transported with minimal twist.
#   6. Camera rolls only around its own forward/centreline axis.
#   7. Roll is deliberately changed near B1/B2/B3...
#   8. Roll is locked between bifurcations.
#   9. Patient remains stationary throughout.
#
# HUD:
#   Fixed clock face.
#   Blue arrow = simulated bronchoscope rotational orientation.
#
# MP4:
#   Uses exactly the SAME applyCameraState() function as playback.
#
# IMPORTANT:
#   This is pre-procedural virtual bronchoscopy.
#   It is NOT real-time bronchoscope localisation.
# =====================================================================

import slicer
import vtk
import qt
import numpy as np

import os
import re
import math
import shutil
import tempfile
import subprocess


# =====================================================================
# SETTINGS
# =====================================================================

ROUTE_NAME = "Route_Centerline"

CAMERA_STEP_MM = 0.50

LOOK_AHEAD_MM = 5.0

VIEW_ANGLE_DEGREES = 85.0

NEAR_CLIPPING_MM = 0.1

FAR_CLIPPING_MM = 1000.0


# =====================================================================
# BIFURCATION BEHAVIOUR
# =====================================================================

ROLL_START_BEFORE_MM = 8.0

ROLL_FINISH_BEFORE_MM = 2.0

PAUSE_BEFORE_MM = 2.0

PAUSE_SECONDS = 1.5

# B1/B2/... markers farther than this from the route are ignored.
MAX_MARKER_DISTANCE_MM = 10.0


# =====================================================================
# PLAYBACK
# =====================================================================

BASE_INTERVAL_MS = 30

DEFAULT_SPEED_PERCENT = 100


# =====================================================================
# VIDEO
# =====================================================================

VIDEO_FPS = 30

VIDEO_CRF = 18

VIDEO_PAUSE_SECONDS = 1.5


# =====================================================================
# ANATOMICAL AXES - SLICER RAS
# =====================================================================

ANATOMICAL_RIGHT = np.array(
    [1.0, 0.0, 0.0],
    dtype=float
)

ANATOMICAL_ANTERIOR = np.array(
    [0.0, 1.0, 0.0],
    dtype=float
)

ANATOMICAL_SUPERIOR = np.array(
    [0.0, 0.0, 1.0],
    dtype=float
)

ANATOMICAL_INFERIOR = np.array(
    [0.0, 0.0, -1.0],
    dtype=float
)


# =====================================================================
# BASIC VECTOR FUNCTIONS
# =====================================================================

def normalize(v):

    v = np.asarray(
        v,
        dtype=float
    )

    n = np.linalg.norm(v)

    if n < 1e-10:
        return None

    return v / n


def rotateAroundAxis(
    vector,
    axis,
    angleRadians
):

    axis = normalize(axis)

    if axis is None:
        return np.asarray(
            vector,
            dtype=float
        ).copy()

    v = np.asarray(
        vector,
        dtype=float
    )

    c = math.cos(
        angleRadians
    )

    s = math.sin(
        angleRadians
    )

    return (
        v * c
        +
        np.cross(
            axis,
            v
        ) * s
        +
        axis
        *
        np.dot(
            axis,
            v
        )
        *
        (1.0 - c)
    )


def smoothstep(x):

    x = max(
        0.0,
        min(
            1.0,
            float(x)
        )
    )

    return (
        x
        *
        x
        *
        (3.0 - 2.0 * x)
    )


def wrapAngle(angle):

    while angle > math.pi:

        angle -= (
            2.0
            *
            math.pi
        )

    while angle < -math.pi:

        angle += (
            2.0
            *
            math.pi
        )

    return angle


def parseClockPosition(description):
    """
    Clock position 1-12 from a marker description.

    Returns None when the description contains no clock value or
    more than one distinct candidate, so that ambiguous text such as
    "branch 6, 3 o'clock" is never silently misread.
    """

    values = {
        int(value)
        for value in re.findall(
            r"\b(1[0-2]|[1-9])\b",
            description or ""
        )
    }

    if len(values) != 1:
        return None

    return values.pop()


# =====================================================================
# FIND ROUTE
# =====================================================================

try:

    routeCurve = slicer.util.getNode(
        ROUTE_NAME
    )

except:

    raise RuntimeError(
        "Route_Centerline was not found."
    )


# =====================================================================
# EXTRACT WORLD COORDINATES
# =====================================================================

polyData = vtk.vtkPolyData()


try:

    routeCurve.GetCurveWorld(
        polyData
    )

except:

    polyData = (
        routeCurve.GetCurveWorld()
    )


if (
    polyData is None
    or
    polyData.GetPoints() is None
):

    raise RuntimeError(
        "Route_Centerline contains no curve geometry."
    )


vtkPoints = (
    polyData.GetPoints()
)


route = np.array(
    [
        vtkPoints.GetPoint(i)

        for i in range(
            vtkPoints.GetNumberOfPoints()
        )
    ],
    dtype=float
)


if len(route) < 2:

    raise RuntimeError(
        "Route_Centerline is too short."
    )


# =====================================================================
# REMOVE DUPLICATE POINTS
# =====================================================================

cleanRoute = [
    route[0]
]


for point in route[1:]:

    if (
        np.linalg.norm(
            point
            -
            cleanRoute[-1]
        )
        >
        0.01
    ):

        cleanRoute.append(
            point
        )


route = np.asarray(
    cleanRoute,
    dtype=float
)


# =====================================================================
# AUTOMATICALLY DETERMINE TRACHEA -> PERIPHERY DIRECTION
#
# The proximal tracheal route should initially travel predominantly
# superior -> inferior.
#
# Compare first and last route portions with anatomical inferior (-Z).
# =====================================================================

testCount = min(
    10,
    len(route) - 1
)


forwardInitialVector = (
    route[testCount]
    -
    route[0]
)


reverseInitialVector = (
    route[-1 - testCount]
    -
    route[-1]
)


forwardInitialVector = normalize(
    forwardInitialVector
)


reverseInitialVector = normalize(
    reverseInitialVector
)


forwardScore = -999.0

reverseScore = -999.0


if forwardInitialVector is not None:

    forwardScore = float(
        np.dot(
            forwardInitialVector,
            ANATOMICAL_INFERIOR
        )
    )


if reverseInitialVector is not None:

    reverseScore = float(
        np.dot(
            reverseInitialVector,
            ANATOMICAL_INFERIOR
        )
    )


routeWasReversed = False


# Routes created by script 03 record their direction explicitly.
# The heuristic is only used for routes without that attribute.
storedDirection = routeCurve.GetAttribute(
    "BBT.RouteDirection"
)

directionIsStored = (
    storedDirection == "StartToTarget"
)


if (
    not directionIsStored
    and
    reverseScore > forwardScore
):

    route = route[::-1].copy()

    routeWasReversed = True


print("")
print(
    "Route direction:",
    (
        "kept as stored (BBT.RouteDirection = StartToTarget)"
        if directionIsStored
        else
        "REVERSED to favour head -> feet"
        if routeWasReversed
        else
        "kept as stored"
    )
)


# =====================================================================
# ROUTE DISTANCE
# =====================================================================

segmentLengths = np.linalg.norm(
    np.diff(
        route,
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


routeLength = float(
    cumulative[-1]
)


if routeLength <= 0:

    raise RuntimeError(
        "Route has zero length."
    )


# =====================================================================
# UNIFORM ROUTE RESAMPLING
# =====================================================================

sampleDistances = np.arange(
    0.0,
    routeLength,
    CAMERA_STEP_MM
)


if (
    len(sampleDistances) == 0
    or
    sampleDistances[-1] < routeLength
):

    sampleDistances = np.append(
        sampleDistances,
        routeLength
    )


flyPoints = []


for distance in sampleDistances:

    index = np.searchsorted(
        cumulative,
        distance
    )

    if index <= 0:

        flyPoints.append(
            route[0]
        )

        continue

    if index >= len(route):

        flyPoints.append(
            route[-1]
        )

        continue


    d0 = cumulative[
        index - 1
    ]

    d1 = cumulative[
        index
    ]


    if abs(
        d1 - d0
    ) < 1e-10:

        fraction = 0.0

    else:

        fraction = (
            distance - d0
        ) / (
            d1 - d0
        )


    point = (
        route[index - 1]
        *
        (1.0 - fraction)
        +
        route[index]
        *
        fraction
    )


    flyPoints.append(
        point
    )


flyPoints = np.asarray(
    flyPoints,
    dtype=float
)


N = len(
    flyPoints
)


if N < 2:

    raise RuntimeError(
        "Insufficient fly-through points."
    )


# =====================================================================
# DISTANCE ALONG FLY ROUTE
# =====================================================================

flySegmentLengths = np.linalg.norm(
    np.diff(
        flyPoints,
        axis=0
    ),
    axis=1
)


flyDistance = np.concatenate(
    (
        [0.0],
        np.cumsum(
            flySegmentLengths
        )
    )
)


# =====================================================================
# SMOOTH ROUTE TANGENTS
#
# These control camera DIRECTION only.
# They do not rotate the patient.
# =====================================================================

tangents = np.zeros_like(
    flyPoints
)


tangentHalfWindow = max(
    1,
    int(
        round(
            LOOK_AHEAD_MM
            /
            CAMERA_STEP_MM
            /
            2.0
        )
    )
)


for i in range(N):

    a = max(
        0,
        i - tangentHalfWindow
    )

    b = min(
        N - 1,
        i + tangentHalfWindow
    )


    if a == b:

        if i < N - 1:

            vector = (
                flyPoints[i + 1]
                -
                flyPoints[i]
            )

        else:

            vector = (
                flyPoints[i]
                -
                flyPoints[i - 1]
            )

    else:

        vector = (
            flyPoints[b]
            -
            flyPoints[a]
        )


    tangent = normalize(
        vector
    )


    if tangent is None:

        if i > 0:

            tangent = (
                tangents[i - 1]
            )

        else:

            tangent = (
                ANATOMICAL_INFERIOR.copy()
            )


    tangents[i] = (
        tangent
    )


# =====================================================================
# INITIAL PATIENT-FIXED CAMERA FRAME
#
# Camera is at proximal trachea.
# Forward follows actual route.
#
# We use PATIENT ANTERIOR (+Y) to establish the initial "up" direction.
#
# This means patient orientation comes from anatomical RAS coordinates,
# not from the current screen orientation.
# =====================================================================

initialForward = (
    tangents[0]
)


initialUp = (
    ANATOMICAL_ANTERIOR
    -
    np.dot(
        ANATOMICAL_ANTERIOR,
        initialForward
    )
    *
    initialForward
)


initialUp = normalize(
    initialUp
)


# If route is almost exactly anterior-posterior,
# use anatomical superior as fallback.
if initialUp is None:

    initialUp = (
        ANATOMICAL_SUPERIOR
        -
        np.dot(
            ANATOMICAL_SUPERIOR,
            initialForward
        )
        *
        initialForward
    )

    initialUp = normalize(
        initialUp
    )


# Final fallback.
if initialUp is None:

    initialUp = (
        ANATOMICAL_RIGHT
        -
        np.dot(
            ANATOMICAL_RIGHT,
            initialForward
        )
        *
        initialForward
    )

    initialUp = normalize(
        initialUp
    )


if initialUp is None:

    raise RuntimeError(
        "Unable to construct initial anatomical camera frame."
    )


# =====================================================================
# MINIMAL-TWIST PARALLEL TRANSPORT
#
# This is the key behaviour:
#
# Patient remains fixed.
#
# Camera forward direction bends with airway.
#
# Camera reference frame is transported with minimum possible twist.
#
# Therefore ordinary airway curvature does NOT create arbitrary roll.
# =====================================================================

transportedUp = np.zeros_like(
    flyPoints
)


transportedUp[0] = (
    initialUp
)


for i in range(
    1,
    N
):

    previousTangent = (
        tangents[i - 1]
    )

    currentTangent = (
        tangents[i]
    )

    previousUp = (
        transportedUp[i - 1]
    )


    rotationAxis = np.cross(
        previousTangent,
        currentTangent
    )


    axisNorm = np.linalg.norm(
        rotationAxis
    )


    tangentDot = float(
        np.clip(
            np.dot(
                previousTangent,
                currentTangent
            ),
            -1.0,
            1.0
        )
    )


    if axisNorm < 1e-9:

        newUp = (
            previousUp.copy()
        )

    else:

        rotationAxis = (
            rotationAxis
            /
            axisNorm
        )


        bendAngle = math.atan2(
            axisNorm,
            tangentDot
        )


        newUp = rotateAroundAxis(
            previousUp,
            rotationAxis,
            bendAngle
        )


    # Remove any numerical component along forward vector.
    newUp = (
        newUp
        -
        np.dot(
            newUp,
            currentTangent
        )
        *
        currentTangent
    )


    newUp = normalize(
        newUp
    )


    if newUp is None:

        newUp = (
            previousUp.copy()
        )


    transportedUp[i] = (
        newUp
    )


# =====================================================================
# FIND CLINICIAN BIFURCATION MARKERS
#
# B1, B2, B3...
#
# Description contains:
#
# 12
# 1
# 2
# ...
# =====================================================================

bifurcations = []


markupNodes = []


for node in slicer.util.getNodesByClass(
    "vtkMRMLMarkupsFiducialNode"
):

    markupNodes.append(
        node
    )


for node in markupNodes:

    nodeName = str(
        node.GetName()
    )


    match = re.fullmatch(
        r"B(\d+)",
        nodeName,
        flags=re.IGNORECASE
    )


    if match is None:

        continue


    if (
        node.GetNumberOfControlPoints()
        <
        1
    ):

        continue


    bifNumber = int(
        match.group(1)
    )


    position = [
        0.0,
        0.0,
        0.0
    ]


    node.GetNthControlPointPositionWorld(
        0,
        position
    )


    position = np.asarray(
        position,
        dtype=float
    )


    description = ""


    try:

        description = str(
            node.GetNthControlPointDescription(
                0
            )
        )

    except:

        pass


    clockValue = parseClockPosition(
        description
    )


    if clockValue is None:

        print(
            "WARNING:",
            nodeName,
            "ignored because its description does not contain "
            "exactly one clock value 1-12."
        )

        continue


    distances = np.linalg.norm(
        flyPoints
        -
        position,
        axis=1
    )


    routeIndex = int(
        np.argmin(
            distances
        )
    )


    markerDistance = float(
        distances[
            routeIndex
        ]
    )


    if markerDistance > MAX_MARKER_DISTANCE_MM:

        print(
            "WARNING:",
            nodeName,
            f"ignored because it is {markerDistance:.1f} mm "
            f"from the route (limit {MAX_MARKER_DISTANCE_MM:.1f} mm)."
        )

        continue


    # ================================================================
    # CLOCK -> CAMERA ROLL
    #
    # 12 = 0°
    # 1  = 30°
    # 2  = 60°
    # ...
    #
    # Camera rotates around its OWN optical axis.
    #
    # Patient does NOT rotate.
    # ================================================================

    clockDegrees = (
        clockValue % 12
    ) * 30.0


    desiredRoll = math.radians(
        -clockDegrees
    )


    bifurcations.append(
        {
            "name": nodeName,

            "number": bifNumber,

            "clock": clockValue,

            "index": routeIndex,

            "distance": float(
                flyDistance[
                    routeIndex
                ]
            ),

            "distanceFromRoute": (
                markerDistance
            ),

            "desiredRoll": (
                desiredRoll
            )
        }
    )


# Anatomical sequence follows route.
bifurcations.sort(
    key=lambda event:
    event["distance"]
)


# =====================================================================
# PRINT BIFURCATION AUDIT
# =====================================================================

print("")
print(
    "=============================================="
)

print(
    "BBT BIFURCATION EVENTS"
)

print(
    "=============================================="
)


if len(
    bifurcations
) == 0:

    print(
        "No valid B1/B2/B3... markers found."
    )

else:

    for event in bifurcations:

        print(
            event["name"],
            "| clock:",
            event["clock"],
            "| route:",
            f"{event['distance']:.1f} mm",
            "| marker-to-route:",
            f"{event['distanceFromRoute']:.1f} mm"
        )


# =====================================================================
# BUILD CAMERA ROLL PROFILE
#
# Roll is zero initially.
#
# Approaching a bifurcation:
#     deliberate roll.
#
# Between bifurcations:
#     roll value remains exactly locked.
# =====================================================================

rollProfile = np.zeros(
    N,
    dtype=float
)


currentLockedRoll = 0.0


for event in bifurcations:

    bifDistance = (
        event["distance"]
    )


    desiredRoll = (
        event["desiredRoll"]
    )


    # Shortest rotational transition.
    deltaRoll = wrapAngle(
        desiredRoll
        -
        currentLockedRoll
    )


    targetRoll = (
        currentLockedRoll
        +
        deltaRoll
    )


    rollStartDistance = max(
        0.0,
        bifDistance
        -
        ROLL_START_BEFORE_MM
    )


    rollFinishDistance = max(
        rollStartDistance,
        bifDistance
        -
        ROLL_FINISH_BEFORE_MM
    )


    startIndex = int(
        np.argmin(
            np.abs(
                flyDistance
                -
                rollStartDistance
            )
        )
    )


    finishIndex = int(
        np.argmin(
            np.abs(
                flyDistance
                -
                rollFinishDistance
            )
        )
    )


    finishIndex = max(
        startIndex,
        finishIndex
    )


    if finishIndex > startIndex:

        denominator = (
            flyDistance[
                finishIndex
            ]
            -
            flyDistance[
                startIndex
            ]
        )


        for i in range(
            startIndex,
            finishIndex + 1
        ):

            if denominator <= 0:

                fraction = 1.0

            else:

                fraction = (
                    flyDistance[i]
                    -
                    flyDistance[
                        startIndex
                    ]
                ) / denominator


            fraction = smoothstep(
                fraction
            )


            rollProfile[i] = (
                currentLockedRoll
                +
                deltaRoll
                *
                fraction
            )


    rollProfile[
        finishIndex:
    ] = targetRoll


    currentLockedRoll = (
        targetRoll
    )


# =====================================================================
# PAUSE EVENTS
# =====================================================================

pauseEvents = []


for event in bifurcations:

    pauseDistance = max(
        0.0,
        event["distance"]
        -
        PAUSE_BEFORE_MM
    )


    pauseIndex = int(
        np.argmin(
            np.abs(
                flyDistance
                -
                pauseDistance
            )
        )
    )


    eventCopy = dict(
        event
    )


    eventCopy[
        "pauseIndex"
    ] = pauseIndex


    pauseEvents.append(
        eventCopy
    )


pauseByIndex = {
    event["pauseIndex"]:
    event

    for event in pauseEvents
}


# =====================================================================
# ACTIVE 3D VIEW
# =====================================================================

layoutManager = (
    slicer.app.layoutManager()
)


if (
    layoutManager is None
    or
    layoutManager.threeDViewCount < 1
):

    raise RuntimeError(
        "No Slicer 3D view available."
    )


threeDWidget = (
    layoutManager.threeDWidget(0)
)


threeDView = (
    threeDWidget.threeDView()
)


viewNode = (
    threeDWidget.mrmlViewNode()
)


cameraNode = (
    slicer.modules.cameras.logic()
    .GetViewActiveCameraNode(
        viewNode
    )
)


if cameraNode is None:

    raise RuntimeError(
        "Could not obtain active camera."
    )


camera = (
    cameraNode.GetCamera()
)


camera.SetViewAngle(
    VIEW_ANGLE_DEGREES
)


# =====================================================================
# IMPORTANT
#
# We manipulate CAMERA ONLY.
#
# No segmentation transform.
# No model transform.
# No volume transform.
# No 3D scene rotation.
#
# Therefore patient remains fixed.
# =====================================================================


# =====================================================================
# HUD RENDERER
# =====================================================================

renderWindow = (
    threeDView.renderWindow()
)


renderer = (
    renderWindow
    .GetRenderers()
    .GetFirstRenderer()
)


# =====================================================================
# REMOVE PREVIOUS BBT HUD IF PRESENT
# =====================================================================

if hasattr(
    slicer,
    "_bbtV7HUDActors"
):

    try:

        for oldActor in (
            slicer._bbtV7HUDActors
        ):

            renderer.RemoveActor2D(
                oldActor
            )

    except:

        pass


hudActors = []


# =====================================================================
# HUD CLOCK
# =====================================================================

HUD_CX = 90

HUD_CY = 90

HUD_RADIUS = 52


circlePoints = vtk.vtkPoints()

circleLines = vtk.vtkCellArray()


CLOCK_SEGMENTS = 72


for i in range(
    CLOCK_SEGMENTS + 1
):

    angle = (
        2.0
        *
        math.pi
        *
        i
        /
        CLOCK_SEGMENTS
    )


    x = (
        HUD_CX
        +
        HUD_RADIUS
        *
        math.sin(
            angle
        )
    )


    y = (
        HUD_CY
        +
        HUD_RADIUS
        *
        math.cos(
            angle
        )
    )


    circlePoints.InsertNextPoint(
        x,
        y,
        0
    )


for i in range(
    CLOCK_SEGMENTS
):

    line = vtk.vtkLine()


    line.GetPointIds().SetId(
        0,
        i
    )


    line.GetPointIds().SetId(
        1,
        i + 1
    )


    circleLines.InsertNextCell(
        line
    )


circlePoly = vtk.vtkPolyData()

circlePoly.SetPoints(
    circlePoints
)

circlePoly.SetLines(
    circleLines
)


circleMapper = vtk.vtkPolyDataMapper2D()

circleMapper.SetInputData(
    circlePoly
)


circleActor = vtk.vtkActor2D()

circleActor.SetMapper(
    circleMapper
)


circleActor.GetProperty().SetColor(
    1.0,
    1.0,
    1.0
)


circleActor.GetProperty().SetLineWidth(
    2.0
)


renderer.AddActor2D(
    circleActor
)


hudActors.append(
    circleActor
)


# =====================================================================
# CLOCK NUMBERS
# =====================================================================

for clockNumber in range(
    1,
    13
):

    angle = math.radians(
        clockNumber
        *
        30.0
    )


    radius = (
        HUD_RADIUS
        +
        16
    )


    x = (
        HUD_CX
        +
        radius
        *
        math.sin(
            angle
        )
    )


    y = (
        HUD_CY
        +
        radius
        *
        math.cos(
            angle
        )
    )


    textActor = vtk.vtkTextActor()


    textActor.SetInput(
        str(
            clockNumber
        )
    )


    textActor.SetPosition(
        x - 6,
        y - 7
    )


    textProperty = (
        textActor
        .GetTextProperty()
    )


    textProperty.SetFontSize(
        14
    )


    textProperty.SetColor(
        1.0,
        1.0,
        1.0
    )


    textProperty.SetBold(
        True
    )


    renderer.AddActor2D(
        textActor
    )


    hudActors.append(
        textActor
    )


# =====================================================================
# BLUE BRONCHOSCOPE ORIENTATION ARROW
# =====================================================================

arrowPoints = vtk.vtkPoints()

arrowLines = vtk.vtkCellArray()


for point in [
    (HUD_CX, HUD_CY, 0),
    (HUD_CX, HUD_CY + 40, 0),
    (HUD_CX, HUD_CY + 40, 0),
    (HUD_CX - 8, HUD_CY + 29, 0),
    (HUD_CX, HUD_CY + 40, 0),
    (HUD_CX + 8, HUD_CY + 29, 0)
]:

    arrowPoints.InsertNextPoint(
        *point
    )


for a, b in [
    (0, 1),
    (2, 3),
    (4, 5)
]:

    line = vtk.vtkLine()


    line.GetPointIds().SetId(
        0,
        a
    )


    line.GetPointIds().SetId(
        1,
        b
    )


    arrowLines.InsertNextCell(
        line
    )


arrowPoly = vtk.vtkPolyData()

arrowPoly.SetPoints(
    arrowPoints
)

arrowPoly.SetLines(
    arrowLines
)


arrowMapper = vtk.vtkPolyDataMapper2D()

arrowMapper.SetInputData(
    arrowPoly
)


arrowActor = vtk.vtkActor2D()

arrowActor.SetMapper(
    arrowMapper
)


arrowActor.GetProperty().SetColor(
    0.1,
    0.55,
    1.0
)


arrowActor.GetProperty().SetLineWidth(
    4.0
)


renderer.AddActor2D(
    arrowActor
)


hudActors.append(
    arrowActor
)


# =====================================================================
# HUD STATUS
# =====================================================================

hudStatus = vtk.vtkTextActor()


hudStatus.SetPosition(
    20,
    170
)


hudStatus.GetTextProperty().SetFontSize(
    16
)


hudStatus.GetTextProperty().SetColor(
    1.0,
    1.0,
    1.0
)


hudStatus.GetTextProperty().SetBold(
    True
)


renderer.AddActor2D(
    hudStatus
)


hudActors.append(
    hudStatus
)


# =====================================================================
# ANATOMICAL REFERENCE HUD
# =====================================================================

anatomyStatus = vtk.vtkTextActor()


anatomyStatus.SetInput(
    "PATIENT FIXED  |  SUPINE  |  camera follows airway"
)


anatomyStatus.SetPosition(
    20,
    195
)


anatomyStatus.GetTextProperty().SetFontSize(
    13
)


anatomyStatus.GetTextProperty().SetColor(
    0.85,
    0.85,
    0.85
)


renderer.AddActor2D(
    anatomyStatus
)


hudActors.append(
    anatomyStatus
)


# =====================================================================
# UPDATE BLUE ARROW
# =====================================================================

def updateArrow(
    rollRadians
):

    displayAngle = (
        -rollRadians
    )


    directionX = math.sin(
        displayAngle
    )


    directionY = math.cos(
        displayAngle
    )


    perpendicularX = (
        -directionY
    )


    perpendicularY = (
        directionX
    )


    tipX = (
        HUD_CX
        +
        40.0
        *
        directionX
    )


    tipY = (
        HUD_CY
        +
        40.0
        *
        directionY
    )


    baseX = (
        HUD_CX
        +
        28.0
        *
        directionX
    )


    baseY = (
        HUD_CY
        +
        28.0
        *
        directionY
    )


    arrowPoints.SetPoint(
        0,
        HUD_CX,
        HUD_CY,
        0
    )


    arrowPoints.SetPoint(
        1,
        tipX,
        tipY,
        0
    )


    arrowPoints.SetPoint(
        2,
        tipX,
        tipY,
        0
    )


    arrowPoints.SetPoint(
        3,
        baseX
        +
        8.0
        *
        perpendicularX,
        baseY
        +
        8.0
        *
        perpendicularY,
        0
    )


    arrowPoints.SetPoint(
        4,
        tipX,
        tipY,
        0
    )


    arrowPoints.SetPoint(
        5,
        baseX
        -
        8.0
        *
        perpendicularX,
        baseY
        -
        8.0
        *
        perpendicularY,
        0
    )


    arrowPoints.Modified()

    arrowPoly.Modified()


# =====================================================================
# SINGLE CAMERA STATE FUNCTION
#
# This is the ONLY function that defines camera pose.
#
# Used by:
#   interactive playback
#   manual slider
#   MP4 rendering
#
# PATIENT GEOMETRY IS NEVER MODIFIED.
# =====================================================================

def applyCameraState(
    index
):

    index = int(
        np.clip(
            index,
            0,
            N - 1
        )
    )


    # ---------------------------------------------------------------
    # CAMERA POSITION
    # ---------------------------------------------------------------

    position = (
        flyPoints[
            index
        ]
    )


    # ---------------------------------------------------------------
    # CAMERA DIRECTION
    #
    # Camera follows actual airway centreline.
    # ---------------------------------------------------------------

    forward = (
        tangents[
            index
        ]
    )


    # ---------------------------------------------------------------
    # MINIMAL-TWIST PATIENT-REFERENCED FRAME
    # ---------------------------------------------------------------

    baseUp = (
        transportedUp[
            index
        ]
    )


    # ---------------------------------------------------------------
    # DELIBERATE BRONCHOSCOPE ROLL
    #
    # Camera rolls.
    # Patient does not.
    # ---------------------------------------------------------------

    roll = float(
        rollProfile[
            index
        ]
    )


    rolledUp = rotateAroundAxis(
        baseUp,
        forward,
        roll
    )


    rolledUp = (
        rolledUp
        -
        np.dot(
            rolledUp,
            forward
        )
        *
        forward
    )


    rolledUp = normalize(
        rolledUp
    )


    if rolledUp is None:

        rolledUp = (
            baseUp
        )


    focalPoint = (
        position
        +
        forward
        *
        LOOK_AHEAD_MM
    )


    # ---------------------------------------------------------------
    # MODIFY CAMERA ONLY
    # ---------------------------------------------------------------

    camera.SetPosition(
        float(
            position[0]
        ),
        float(
            position[1]
        ),
        float(
            position[2]
        )
    )


    camera.SetFocalPoint(
        float(
            focalPoint[0]
        ),
        float(
            focalPoint[1]
        ),
        float(
            focalPoint[2]
        )
    )


    camera.SetViewUp(
        float(
            rolledUp[0]
        ),
        float(
            rolledUp[1]
        ),
        float(
            rolledUp[2]
        )
    )


    camera.OrthogonalizeViewUp()


    camera.SetViewAngle(
        VIEW_ANGLE_DEGREES
    )


    camera.SetClippingRange(
        NEAR_CLIPPING_MM,
        FAR_CLIPPING_MM
    )


    cameraNode.Modified()


    updateArrow(
        roll
    )


    # ---------------------------------------------------------------
    # BIFURCATION HUD
    # ---------------------------------------------------------------

    eventText = ""


    for event in bifurcations:

        distanceAhead = (
            event["distance"]
            -
            flyDistance[
                index
            ]
        )


        if (
            -2.0
            <=
            distanceAhead
            <=
            ROLL_START_BEFORE_MM + 2.0
        ):

            eventText = (
                event["name"]
                +
                "  ->  "
                +
                str(
                    event["clock"]
                )
                +
                " o'clock"
            )

            break


    hudStatus.SetInput(
        eventText
    )


    threeDView.forceRender()


# =====================================================================
# INITIAL STATE
# =====================================================================

currentIndex = 0

playing = False

waitingAtBifurcation = False

pauseAlreadyTriggered = set()

recording = False

cancelRecording = False


applyCameraState(
    0
)


# =====================================================================
# REMOVE OLD CONTROLLER
# =====================================================================

if hasattr(
    slicer,
    "_bbtV7Controller"
):

    try:

        slicer._bbtV7Controller.close()

    except:

        pass


if hasattr(
    slicer,
    "_bbtV7Timer"
):

    try:

        slicer._bbtV7Timer.stop()

    except:

        pass


# =====================================================================
# CONTROLLER WINDOW
# =====================================================================

controller = qt.QWidget()


controller.setWindowTitle(
    "BBT Virtual Bronchoscopy v7"
)


controller.setWindowFlags(
    qt.Qt.Window
)


controller.setMinimumWidth(
    580
)


controller.setEnabled(
    True
)


layout = qt.QVBoxLayout(
    controller
)


title = qt.QLabel(
    "<b>BBT Virtual Bronchoscopy — Patient-Fixed Camera</b>"
)


layout.addWidget(
    title
)


orientationLabel = qt.QLabel(
    "Patient fixed supine | initial reference: head →feet | "
    "camera rolls around centreline"
)


orientationLabel.wordWrap = True


layout.addWidget(
    orientationLabel
)


infoLabel = qt.QLabel(
    "Route: "
    +
    f"{routeLength:.1f} mm"
    +
    "   |   Bifurcations: "
    +
    str(
        len(
            bifurcations
        )
    )
)


layout.addWidget(
    infoLabel
)


statusLabel = qt.QLabel(
    "Ready at proximal route"
)


layout.addWidget(
    statusLabel
)


# =====================================================================
# CONTROL BUTTONS
# =====================================================================

buttonLayout = qt.QHBoxLayout()


playButton = qt.QPushButton(
    "▶Play"
)


pauseButton = qt.QPushButton(
    "ⅡPause"
)


stopButton = qt.QPushButton(
    "■Stop"
)


restartButton = qt.QPushButton(
    "↻Restart"
)


for button in [
    playButton,
    pauseButton,
    stopButton,
    restartButton
]:

    button.setEnabled(
        True
    )

    buttonLayout.addWidget(
        button
    )


layout.addLayout(
    buttonLayout
)


# =====================================================================
# POSITION SLIDER
# =====================================================================

positionLayout = qt.QHBoxLayout()


positionLabel = qt.QLabel(
    "Position"
)


positionSlider = qt.QSlider(
    qt.Qt.Horizontal
)


positionSlider.minimum = 0

positionSlider.maximum = (
    N - 1
)

positionSlider.value = 0


positionLayout.addWidget(
    positionLabel
)


positionLayout.addWidget(
    positionSlider
)


layout.addLayout(
    positionLayout
)


# =====================================================================
# SPEED
# =====================================================================

speedLayout = qt.QHBoxLayout()


speedLabel = qt.QLabel(
    "Speed"
)


speedSlider = qt.QSlider(
    qt.Qt.Horizontal
)


speedSlider.minimum = 25

speedSlider.maximum = 300

speedSlider.value = (
    DEFAULT_SPEED_PERCENT
)


speedDisplay = qt.QLabel(
    "1.00×"
)


speedLayout.addWidget(
    speedLabel
)


speedLayout.addWidget(
    speedSlider
)


speedLayout.addWidget(
    speedDisplay
)


layout.addLayout(
    speedLayout
)


# =====================================================================
# BIFURCATION SUMMARY
# =====================================================================

bifLabel = qt.QLabel(
    "<b>Bifurcation orientation</b>"
)


layout.addWidget(
    bifLabel
)


bifText = qt.QLabel()

bifText.wordWrap = True


if bifurcations:

    bifText.text = "   |   ".join(
        [
            (
                event["name"]
                +
                ": "
                +
                str(
                    event["clock"]
                )
                +
                " o'clock"
            )

            for event in bifurcations
        ]
    )

else:

    bifText.text = (
        "No B1/B2/B3... markers with clock descriptions found."
    )


layout.addWidget(
    bifText
)


# =====================================================================
# VIDEO CONTROLS
# =====================================================================

videoLayout = qt.QHBoxLayout()


saveButton = qt.QPushButton(
    "●Save MP4"
)


cancelButton = qt.QPushButton(
    "Cancel MP4"
)


cancelButton.enabled = False


videoLayout.addWidget(
    saveButton
)


videoLayout.addWidget(
    cancelButton
)


layout.addLayout(
    videoLayout
)


# =====================================================================
# TIMERS
# =====================================================================

timer = qt.QTimer()

pauseTimer = qt.QTimer()


pauseTimer.setSingleShot(
    True
)


def updateTimerSpeed():

    speed = (
        float(
            speedSlider.value
        )
        /
        100.0
    )


    speedDisplay.text = (
        f"{speed:.2f}×"
    )


    interval = max(
        1,
        int(
            round(
                BASE_INTERVAL_MS
                /
                speed
            )
        )
    )


    timer.setInterval(
        interval
    )


updateTimerSpeed()


# =====================================================================
# RESUME AFTER BIFURCATION
# =====================================================================

def resumeAfterBifurcation():

    global waitingAtBifurcation


    waitingAtBifurcation = False


    if playing:

        timer.start()

        statusLabel.text = (
            "Playing"
        )


pauseTimer.connect(
    "timeout()",
    resumeAfterBifurcation
)


# =====================================================================
# PLAYBACK STEP
# =====================================================================

def timerStep():

    global currentIndex
    global playing
    global waitingAtBifurcation


    if not playing:

        return


    if waitingAtBifurcation:

        return


    nextIndex = (
        currentIndex + 1
    )


    if nextIndex >= N:

        playing = False

        timer.stop()

        statusLabel.text = (
            "Route target reached"
        )

        return


    currentIndex = (
        nextIndex
    )


    applyCameraState(
        currentIndex
    )


    positionSlider.blockSignals(
        True
    )


    positionSlider.value = (
        currentIndex
    )


    positionSlider.blockSignals(
        False
    )


    # ---------------------------------------------------------------
    # AUTOMATIC PAUSE
    # ---------------------------------------------------------------

    if (
        currentIndex in pauseByIndex
        and
        currentIndex
        not in pauseAlreadyTriggered
    ):

        event = (
            pauseByIndex[
                currentIndex
            ]
        )


        pauseAlreadyTriggered.add(
            currentIndex
        )


        waitingAtBifurcation = True


        timer.stop()


        statusLabel.text = (
            "Paused before "
            +
            event["name"]
            +
            " — "
            +
            str(
                event["clock"]
            )
            +
            " o'clock"
        )


        pauseTimer.start(
            int(
                PAUSE_SECONDS
                *
                1000
            )
        )


        return


    percentage = (
        100.0
        *
        currentIndex
        /
        max(
            1,
            N - 1
        )
    )


    statusLabel.text = (
        "Playing — "
        +
        f"{percentage:.0f}%"
    )


timer.connect(
    "timeout()",
    timerStep
)


# =====================================================================
# PLAY
# =====================================================================

def play():

    global playing


    playing = True


    updateTimerSpeed()


    timer.start()


    statusLabel.text = (
        "Playing"
    )


# =====================================================================
# PAUSE
# =====================================================================

def pause():

    global playing


    playing = False


    timer.stop()

    pauseTimer.stop()


    statusLabel.text = (
        "Paused"
    )


# =====================================================================
# RESET
# =====================================================================

def resetToStart():

    global currentIndex
    global waitingAtBifurcation
    global pauseAlreadyTriggered


    timer.stop()

    pauseTimer.stop()


    currentIndex = 0

    waitingAtBifurcation = False

    pauseAlreadyTriggered = set()


    positionSlider.blockSignals(
        True
    )


    positionSlider.value = 0


    positionSlider.blockSignals(
        False
    )


    applyCameraState(
        0
    )


# =====================================================================
# STOP
# =====================================================================

def stop():

    global playing


    playing = False


    resetToStart()


    statusLabel.text = (
        "Stopped — proximal route"
    )


# =====================================================================
# RESTART
# =====================================================================

def restart():

    global playing


    playing = False


    resetToStart()


    play()


# =====================================================================
# POSITION SLIDER
# =====================================================================

def positionChanged(
    value
):

    global currentIndex
    global playing


    if recording:

        return


    playing = False


    timer.stop()

    pauseTimer.stop()


    currentIndex = int(
        value
    )


    applyCameraState(
        currentIndex
    )


    percentage = (
        100.0
        *
        currentIndex
        /
        max(
            1,
            N - 1
        )
    )


    statusLabel.text = (
        "Manual position — "
        +
        f"{percentage:.0f}%"
    )


# =====================================================================
# SPEED
# =====================================================================

def speedChanged(
    value
):

    updateTimerSpeed()


# =====================================================================
# FIND FFMPEG
# =====================================================================

def findFFmpeg():

    candidates = []


    detected = shutil.which(
        "ffmpeg"
    )


    if detected:

        candidates.append(
            detected
        )


    candidates.extend(
        [
            "/opt/homebrew/bin/ffmpeg",
            "/usr/local/bin/ffmpeg",
            "/opt/local/bin/ffmpeg",
            "/usr/bin/ffmpeg"
        ]
    )


    for candidate in candidates:

        if (
            candidate
            and
            os.path.isfile(
                candidate
            )
            and
            os.access(
                candidate,
                os.X_OK
            )
        ):

            return candidate


    return None


# =====================================================================
# CAPTURE FRAME
# =====================================================================

def captureFrame(
    filename
):

    renderWindow.Render()


    capture = vtk.vtkWindowToImageFilter()


    capture.SetInput(
        renderWindow
    )


    capture.SetInputBufferTypeToRGB()


    capture.ReadFrontBufferOff()


    capture.Update()


    writer = vtk.vtkPNGWriter()


    writer.SetFileName(
        filename
    )


    writer.SetInputConnection(
        capture.GetOutputPort()
    )


    writer.Write()


# =====================================================================
# CANCEL MP4
# =====================================================================

def cancelMP4():

    global cancelRecording


    if recording:

        cancelRecording = True


        statusLabel.text = (
            "Cancelling MP4..."
        )


# =====================================================================
# SAVE MP4
#
# applyCameraState() is used directly.
#
# No separate camera algorithm exists for export.
# =====================================================================

def saveMP4():

    global recording
    global cancelRecording
    global playing


    playing = False


    timer.stop()

    pauseTimer.stop()


    ffmpeg = findFFmpeg()


    if ffmpeg is None:

        qt.QMessageBox.critical(
            controller,
            "FFmpeg not found",
            (
                "FFmpeg was not found.\n\n"
                "Install FFmpeg and try again."
            )
        )

        return


    defaultPath = os.path.expanduser(
        "~/Desktop/BBT_VirtualBronchoscopy_v7.mp4"
    )


    outputPath = qt.QFileDialog.getSaveFileName(
        controller,
        "Save virtual bronchoscopy",
        defaultPath,
        "MP4 video (*.mp4)"
    )


    if isinstance(
        outputPath,
        (tuple, list)
    ):

        outputPath = (
            outputPath[0]
            if len(outputPath) > 0
            else ""
        )


    outputPath = str(
        outputPath
    )


    if not outputPath:

        return


    if not outputPath.lower().endswith(
        ".mp4"
    ):

        outputPath += (
            ".mp4"
        )


    tempDirectory = tempfile.mkdtemp(
        prefix="BBT_v7_"
    )


    recording = True

    cancelRecording = False


    saveButton.enabled = False

    cancelButton.enabled = True


    playButton.enabled = False

    pauseButton.enabled = False

    stopButton.enabled = False

    restartButton.enabled = False

    positionSlider.enabled = False


    try:

        frameNumber = 0


        for index in range(N):

            if cancelRecording:

                raise RuntimeError(
                    "MP4 export cancelled."
                )


            # =========================================================
            # SAME CAMERA STATE AS INTERACTIVE NAVIGATION
            # =========================================================

            applyCameraState(
                index
            )


            slicer.app.processEvents()


            framePath = os.path.join(
                tempDirectory,
                f"frame_{frameNumber:06d}.png"
            )


            captureFrame(
                framePath
            )


            frameNumber += 1


            # =========================================================
            # SAME BIFURCATION PAUSE
            # =========================================================

            if index in pauseByIndex:

                event = (
                    pauseByIndex[
                        index
                    ]
                )


                duplicateFrames = int(
                    round(
                        VIDEO_PAUSE_SECONDS
                        *
                        VIDEO_FPS
                    )
                )


                statusLabel.text = (
                    "MP4 — pause before "
                    +
                    event["name"]
                )


                for dummy in range(
                    duplicateFrames
                ):

                    if cancelRecording:

                        raise RuntimeError(
                            "MP4 export cancelled."
                        )


                    framePath = os.path.join(
                        tempDirectory,
                        f"frame_{frameNumber:06d}.png"
                    )


                    captureFrame(
                        framePath
                    )


                    frameNumber += 1


            percentage = (
                100.0
                *
                index
                /
                max(
                    1,
                    N - 1
                )
            )


            statusLabel.text = (
                "Rendering MP4 — "
                +
                f"{percentage:.0f}%"
            )


            slicer.app.processEvents()


        # =============================================================
        # FFMPEG
        # =============================================================

        statusLabel.text = (
            "Encoding MP4..."
        )


        slicer.app.processEvents()


        inputPattern = os.path.join(
            tempDirectory,
            "frame_%06d.png"
        )


        command = [
            ffmpeg,

            "-y",

            "-framerate",
            str(
                VIDEO_FPS
            ),

            "-i",
            inputPattern,

            "-vf",
            "scale=trunc(iw/2)*2:trunc(ih/2)*2",

            "-c:v",
            "libx264",

            "-preset",
            "medium",

            "-crf",
            str(
                VIDEO_CRF
            ),

            "-pix_fmt",
            "yuv420p",

            "-movflags",
            "+faststart",

            outputPath
        ]


        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )


        if result.returncode != 0:

            print(
                result.stderr
            )


            raise RuntimeError(
                "FFmpeg encoding failed."
            )


        statusLabel.text = (
            "MP4 saved"
        )


        qt.QMessageBox.information(
            controller,
            "BBT",
            (
                "MP4 successfully saved:\n\n"
                +
                outputPath
            )
        )


    except Exception as error:

        statusLabel.text = (
            "MP4 export stopped"
        )


        qt.QMessageBox.warning(
            controller,
            "BBT MP4",
            str(
                error
            )
        )


    finally:

        try:

            shutil.rmtree(
                tempDirectory
            )

        except:

            pass


        recording = False

        cancelRecording = False


        saveButton.enabled = True

        cancelButton.enabled = False


        playButton.enabled = True

        pauseButton.enabled = True

        stopButton.enabled = True

        restartButton.enabled = True

        positionSlider.enabled = True


# =====================================================================
# CONNECT CONTROLS
# =====================================================================

playButton.connect(
    "clicked()",
    play
)


pauseButton.connect(
    "clicked()",
    pause
)


stopButton.connect(
    "clicked()",
    stop
)


restartButton.connect(
    "clicked()",
    restart
)


positionSlider.connect(
    "valueChanged(int)",
    positionChanged
)


speedSlider.connect(
    "valueChanged(int)",
    speedChanged
)


saveButton.connect(
    "clicked()",
    saveMP4
)


cancelButton.connect(
    "clicked()",
    cancelMP4
)


# =====================================================================
# KEEP PYTHONQT OBJECTS ALIVE
# =====================================================================

slicer._bbtV7Controller = (
    controller
)

slicer._bbtV7Timer = (
    timer
)

slicer._bbtV7PauseTimer = (
    pauseTimer
)

slicer._bbtV7HUDActors = (
    hudActors
)

slicer._bbtV7FlyPoints = (
    flyPoints
)

slicer._bbtV7Tangents = (
    tangents
)

slicer._bbtV7TransportedUp = (
    transportedUp
)

slicer._bbtV7RollProfile = (
    rollProfile
)


# =====================================================================
# SHOW WINDOW
# =====================================================================

controller.show()

controller.raise_()

controller.activateWindow()


applyCameraState(
    0
)


# =====================================================================
# FINAL REPORT
# =====================================================================

print("")
print(
    "=============================================="
)

print(
    "BBT VIRTUAL BRONCHOSCOPY v7 READY"
)

print(
    "=============================================="
)

print(
    f"Route length: {routeLength:.1f} mm"
)

print(
    f"Camera positions: {N}"
)

print(
    f"Bifurcations: {len(bifurcations)}"
)

print(
    "Patient frame: fixed Slicer RAS"
)

print(
    "Patient position: supine reference"
)

print(
    "Initial travel preference: head -> feet"
)

print(
    "Camera: centreline-following"
)

print(
    "Frame: minimal-twist parallel transport"
)

print(
    "Roll: camera only"
)

print(
    "Patient/model rotation: NONE"
)

print(
    "Clock HUD: fixed"
)

print(
    "Blue arrow: bronchoscope orientation"
)

print(
    "MP4 camera model: identical to interactive playback"
)

print("")
