# CT-based virtual bronchoscopy and bronchial branch tracing workflow

**Step-by-step README for reproducing the current research prototype**\
**Version 1.1 \| Research prototype \| 6 October 2026**

## Purpose

This README explains how to reproduce the CT-based airway
reconstruction, route creation and virtual bronchoscopy workflow
developed so far in 3D Slicer.

The workflow takes a chest CT, marks a pulmonary lesion, segments the
visible airway tree, defines a route from the trachea to the segmented
airway nearest the lesion, displays the route in 3D, and creates a
virtual bronchoscopy fly-through.

This is a **research and pre-procedural planning prototype**. It is not
a medical device and it does not provide live bronchoscope tracking. The
virtual bronchoscopy shows a path through the CT-derived airway model.
It does **not** prove that the real bronchoscope, radial EBUS probe or
biopsy tool is in the same position during a procedure.

### Anatomical rule

The CT is the anatomical source of truth.

Do not: - invent a distal airway that cannot be seen on the CT; -
connect two separated airway segments merely to obtain a continuous
route; - change airway diameter, direction or branch geometry to make
navigation easier; - accept a smooth-looking 3D model as evidence that
an airway exists.

Any manual correction must be supported by the CT and recorded.

------------------------------------------------------------------------

## 1. Files required

Before starting, assemble one folder containing:

1.  The source chest CT in DICOM format.
2.  The Python scripts used by this workflow:
    -   `1 BBT CT Series Selector v1.0`
    -   `2 Make 3D visible`
    -   `3 BBT Continuous Smooth Route v4.0`
    -   `4 BBT Virtual Bronch v7.0`
3.  `BBT_Slicer_SKILL.md`, if you are using an AI assistant to help with
    Slicer or Python.
4.  A saved demonstration Slicer scene, if available.
5.  A short record of the software versions used.
6.  A known expected result for the demonstration case, including the
    intended branch sequence, final CT-visible airway and expected
    output files.

For research transfer, use a de-identified CT for which you have
permission to share.

**Reproduction check:** another operator should be able to open the
demonstration case and reproduce the same CT-supported branch sequence
and fly-through. A visually similar 3D model is not enough if the
branches differ.

------------------------------------------------------------------------

## 2. Install and prepare 3D Slicer

The final verified workflow used **3D Slicer 5.12.4, revision 34654**,
on a 2024 iMac M4 running macOS Tahoe 26.6.2. Other recent versions may
work, but they should be recorded because menus, extensions and Python
interfaces can change.

### 2.1 Learn the few parts of Slicer you need

You do not need to learn the whole program. Before starting, locate
these controls:

-   **DICOM**: imports CT examinations.
-   **Data**: lists all loaded CT volumes, markups, segmentations and
    models.
-   **Modules** and its search box: opens Slicer tools.
-   **Markups**: creates points and regions of interest.
-   **Save**: saves the Slicer scene and outputs.
-   **Extensions Manager**: installs additional modules.
-   **Python Console**: runs the supplied Python scripts.

If the Python Console is not visible, open it from **View \> Python
Interactor**.

### 2.2 Install the extensions

Open **Extensions Manager**, install the required extensions, then
restart Slicer.

Versions used in the verified setup were:

-   SlicerVMTK: `bd2332b` (17 September 2026)
-   PyTorch: `74f8216` (10 September 2026)
-   TotalSegmentator: `270cac2` (10 September 2026)
-   NNUNet: `0cb736d` (10 September 2026)

If you use different versions, write them down. A later extension
version can behave differently from the version used to create this
workflow.

### 2.3 Prepare folders

Create separate folders for:

-   original DICOM;
-   working Slicer scenes;
-   reviewed outputs;
-   exported images and videos.

Do not send identifiable clinical DICOM, CT screenshots or other patient
data to a public generative AI service.

### 2.4 Optional AI coding assistant

An AI chatbot can help explain Slicer errors and draft or revise Python
code. If used, attach `BBT_Slicer_SKILL.md` to the AI conversation so
that the assistant has the project-specific context.

The prototype was developed by repeatedly: 1. describing the required
Slicer behaviour; 2. generating or revising Python; 3. running the code
locally in Slicer; 4. checking the result against the source CT; 5.
correcting the script when the result was wrong.

AI-generated code is **not automatically correct**. Always test it in
Slicer and verify every anatomical output against the CT.

------------------------------------------------------------------------

## 3. Import the CT and select the correct series

### 3.1 Import

1.  Open 3D Slicer.
2.  Open the **DICOM** module.
3.  Import the patient's DICOM folder.
4.  Load the chest CT series that you want to inspect.

### 3.2 Select the best CT series

After importing the available CT series, paste and run:

`1 BBT CT Series Selector v1.0`

in the Slicer Python Console.

The script examines scalar CT volumes loaded in the scene and ranks
likely candidates for airway segmentation. Treat its ranking as an aid,
not as an automatic final decision.

For manual selection, prefer: - a contiguous diagnostic chest CT; - the
thinnest available slices; - a reconstruction that preserves small
airway detail; - a lung or relatively sharp reconstruction kernel,
provided noise does not make the airway lumen difficult to follow.

Thin-section CT is particularly important for bronchial branch tracing.
Published manual branch-mapping methods commonly use thin CT sections,
and performance deteriorates as distal airways become smaller.

### 3.3 Check the CT before proceeding

1.  Display the CT using a lung window preset.
2.  Inspect the **axial, coronal and sagittal** planes.
3.  Confirm that you can identify the trachea, target lobe and lesion.
4.  Scroll through the intended route and look for motion artefact,
    missing slices or severe noise.
5.  Use one chosen CT series as the anatomical reference for the case.

**Stop here** if the CT does not show enough anatomy to support the
intended branch sequence.

------------------------------------------------------------------------

## 4. Mark the lesion

1.  Select the **Red** slice viewer, normally the axial view.
2.  Right-click in the image and choose the **CT-Lung** window/level
    preset.
3.  Enlarge the axial viewer if this makes the lesion easier to see.
4.  Scroll to the lesion.
5.  Open **Markups** and create a **Point List**.
6.  Place one point in the lesion.
7.  Check the point in axial, coronal and sagittal views. Move it until
    its position is correct in all three planes.
8.  Rename the **markup node** `Lesion`.
9.  If the individual control point has a separate label, also label it
    `Lesion`.

The exact node name matters because the supplied Python scripts search
for named Slicer nodes. A spelling or capitalisation mismatch can
produce an error such as `no node was found`.

------------------------------------------------------------------------

## 5. Decide whether to crop the CT

Cropping can reduce processing time, but it can also remove anatomy that
the route algorithm needs.

### Preferred reproducible approach

For the first reproduction attempt, either:

-   run airway segmentation on the **full selected chest CT**, or
-   crop a volume that still contains the **entire intended connected
    route from the trachea to the target lobe and distal target
    airway**.

Do **not** crop only a small box around the lesion if you intend to
calculate a continuous centreline from the trachea. A local crop around
the lesion cannot contain the trachea, so it cannot by itself produce a
trachea-to-target route.

A small local crop may later be used experimentally to try to recover
more distal airway detail, but combining that local result with the main
airway tree requires an explicit, CT-verified merging step. That merging
step is not yet part of this README.

### 5.1 Create a whole-route crop, if required

1.  Open **Markups \> ROI** and create a new region of interest (ROI).
2.  Resize the ROI in axial, coronal and sagittal views.
3.  Make sure it includes:
    -   the tracheal starting point;
    -   the main and lobar bronchi on the intended side;
    -   the complete expected airway route;
    -   the lesion;
    -   enough surrounding tissue to avoid cutting an airway at the ROI
        boundary.
4.  Search for and open **Crop Volume**.
5.  Set:
    -   **Input volume:** selected thin-section CT
    -   **Input ROI:** the ROI just created
    -   **Output volume:** create a new volume, for example
        `CT_ROUTE_ROI`
    -   **Voxel-based cropping:** ON, if available
    -   **Isotropic voxel size:** OFF
    -   **Spacing scale:** 1.0
    -   **Interpolation:** use the normal image-volume interpolation
        option; do not deliberately reduce the CT resolution
6.  Click **Apply**.
7.  Inspect the new cropped volume and confirm that the full route
    remains present.
8.  Hide the ROI box from the Data panel if it obstructs the view.

**Checkpoint:** the trachea, intended route and lesion must all be
visible inside the volume that will be sent to TotalSegmentator.

------------------------------------------------------------------------

## 6. Generate the airway segmentation

1.  Search **Modules** for **TotalSegmentator** and open it.
2.  Set **Input volume** to the selected full CT or the whole-route
    cropped CT.
3.  Select the task that produces the lung airway segmentation. In
    current TotalSegmentator terminology this is the **`lung_vessels`**
    task, which can output lung arteries, lung veins, lung airways and
    airway wall. The wording shown by the Slicer extension may differ
    slightly by version.
4.  Use the normal/full-resolution mode rather than a deliberately
    reduced-resolution fast mode when distal airway detail is the
    priority.
5.  Create a new output segmentation.
6.  If available in your installed version, **Use standard segment
    names** can be enabled.
7.  **Higher-order resampling** may produce smoother upsampling, but it
    should not be treated as proof of better anatomical accuracy.
8.  Click **Apply** and wait for the process to finish.

CPU processing can take several minutes or longer. Runtime depends on CT
size, crop size, computer hardware and software version.

**Do not ignore an error merely because a segmentation appears on
screen.** Read the error. Confirm that the expected airway segment was
produced and that the process completed. If you use an AI assistant to
interpret an error, provide the text of the error, not patient data.

------------------------------------------------------------------------

## 7. Inspect and correct the segmentation

This is the most important quality-control step.

1.  Open **Segment Editor** or **Data**.
2.  Hide pulmonary artery and pulmonary vein segments if they obstruct
    the airway view.
3.  Display the airway segmentation over the CT in axial, coronal and
    sagittal planes.
4.  Show the airway surface in the 3D viewer.
5.  Inspect:
    -   trachea;
    -   main bronchi;
    -   lobar bronchi;
    -   every branch on the intended route;
    -   the most distal CT-supported branch near the lesion.
6.  Look specifically for:
    -   leakage outside the airway lumen;
    -   false connections between neighbouring branches;
    -   missing airway segments;
    -   holes or abrupt breaks;
    -   a distal branch shown in 3D that is not supported by the CT.

### Manual correction rule

Correct a segmentation only when the source CT supports the correction.

Keep the original segmentation or an unedited Slicer scene before manual
editing. Record what you changed.

Do not bridge a gap merely because two pieces appear likely to belong
together. If a gap is not visible on CT, the route must stop or the
uncertainty must be recorded.

------------------------------------------------------------------------

## 8. Mark the start and end of the route

The route script needs two named points.

### RouteTarget

1.  In the 3D view or slice views, create a markup point inside the
    **most distal segmented airway that is supported by CT and is
    appropriate for approaching the lesion**.
2.  Check the point in axial, coronal and sagittal CT.
3.  Make sure the point lies inside the airway segmentation, not simply
    near it.
4.  Name the markup node `RouteTarget`.

`RouteTarget` is an airway endpoint for route calculation. It is not
necessarily the centre of the lesion.

### RouteStart

1.  Place a second markup point inside the tracheal lumen, within the
    connected airway segmentation.
2.  Check its position on CT.
3.  Name the markup node `RouteStart`.

**Checkpoint:** both points must lie inside the same connected,
CT-supported airway tree. If they are in disconnected components, a
continuous route cannot be calculated without introducing an unsupported
connection.

------------------------------------------------------------------------

## 9. Display the 3D airway model

Run:

`2 Make 3D visible`

in the Slicer Python Console.

Confirm that the airway tree appears correctly in the 3D viewer. Use the
**Data** panel to show or hide nodes with the eye icon.

Changing display properties such as point size, label size or colour is
acceptable. Do not change anatomical geometry merely to improve
appearance.

------------------------------------------------------------------------

## 10. Create the centreline route

Run:

`3 BBT Continuous Smooth Route v4.0`

in the Slicer Python Console.

The script should create a route through the segmented airway from
`RouteStart` towards `RouteTarget`.

Check the result before doing anything else:

1.  Confirm that the route begins at `RouteStart`.
2.  Confirm that it remains inside the segmented airway.
3.  Inspect every major turn and bifurcation against the CT.
4.  Confirm that it ends at or near `RouteTarget`.
5.  Look for shortcuts across airway walls or jumps between branches.
6.  If the route is wrong, do not repair it by inventing anatomy.
    Recheck the segmentation, point placement and connectivity.

A centreline is a mathematical path through the segmentation. It is only
as anatomically valid as the segmentation from which it was calculated.

### Built-in safety checks

The start and target points are chosen from a list of all markup nodes;
the names above are a convention. The script then:

-   stops if either point is more than 5 mm from the selected airway
    segment (`MAX_START_SNAP_MM`, `MAX_TARGET_SNAP_MM`), instead of
    snapping it to an unrelated airway voxel;
-   only proposes bridges to fragments that lie closer to the target
    than the start-connected airway already reaches, so isolated
    fragments behind or beside the route are never bridged;
-   centres the slice views on each proposed gap before asking for
    authorisation;
-   constrains the smoothed route to the Start-to-Target connected
    component only, not to every fragment of the navigation mask;
-   reports how many final route points lie outside that component and
    stores the count in `BBT.PointsOutsideRouteComponent`.

Projection keeps points inside the airway but places corrected points
at the wall, not at the lumen centre. A non-zero outside count, or any
projected points near a bifurcation, requires review on CT.

------------------------------------------------------------------------

## 11. Create the virtual bronchoscopy

Run:

`4 BBT Virtual Bronch v7.0`

in the Slicer Python Console.

The current script moves a virtual camera along the calculated route and
is designed to reproduce a bronchoscope-like view, including controlled
camera roll and pauses around branch points.

Review the complete fly-through before export.

At each bifurcation: 1. pause the animation; 2. compare the virtual view
with the 3D model; 3. compare the branch decision with the source CT; 4.
confirm that the selected branch is the intended CT-supported branch.

A virtual bronchoscopy view is a rendering of the CT-derived model. It
is **not** the same as the live endoscopic view and is not real-time
localisation.

### Route direction

Routes created by script 3 carry `BBT.RouteDirection = StartToTarget`
and are played in that order. Other curves named `Route_Centerline`
fall back to a heuristic that prefers initial head-to-feet travel;
check the "Route direction" line in the console.

### Bifurcation markers (optional)

Before running script 4, a point can be placed at each bifurcation:

1.  Create one markup point list per bifurcation and name the **node**
    `B1`, `B2`, `B3`…; events are ordered by their position along the
    route, not by number.
2.  Place the first control point at the bifurcation, within 10 mm of
    the route (`MAX_MARKER_DISTANCE_MM`); farther markers are ignored.
3.  Set the description of that control point to the desired clock
    position only, for example `3`. Descriptions with no value, or with
    more than one different number from 1 to 12, are ignored with a
    warning.

The camera roll convention is described in
`scripts/BBT_Slicer_SKILL.md`. Markers are teaching events chosen by the
operator; they are not automatic anatomical labels.

The scripts do not read the lesion markup; it is a reference for the
operator when placing `RouteTarget`.

------------------------------------------------------------------------

## 12. Export and save

Save the Slicer scene before exporting video.

The review package should contain, where available:

-   the saved Slicer scene;
-   the airway segmentation;
-   the route/centreline;
-   the lesion, `RouteStart` and `RouteTarget` markups;
-   representative CT screenshots;
-   representative bifurcation images;
-   the virtual bronchoscopy video;
-   a record of software and script versions;
-   a record of manual corrections and unresolved uncertainty.

The 3D reconstruction or fly-through can be exported as an MP4 if the
recording/export part of the script and local video dependencies are
working. If MP4 export fails, preserve the Slicer scene and still
images, fix the export problem, and repeat the recording. Do not treat
an incomplete or corrupted video as a final output.

------------------------------------------------------------------------

## 13. Failure handling

  -----------------------------------------------------------------------
  Problem                             What to do
  ----------------------------------- -----------------------------------
  No airway mask or an empty mask     Check the CT input,
                                      TotalSegmentator task, model
                                      download, extension version and
                                      Slicer error log. Do not proceed.

  Trachea is missing after cropping   The ROI is too small for a
                                      trachea-to-target route. Return to
                                      the full CT or create a larger
                                      whole-route ROI.

  Distal airway is missing            Recheck the CT in all three planes.
                                      If the airway is visible but missed
                                      by segmentation, perform a
                                      documented manual correction or a
                                      separately validated local
                                      re-segmentation. If it is not
                                      visible on CT, do not invent it.

  False distal connection             Inspect the CT in all planes and
                                      remove or reject the false
                                      connection.

  `RouteStart` and `RouteTarget` are  Recheck segmentation and point
  disconnected                        placement. Do not create an
                                      invisible connecting airway merely
                                      to make the route run.

  Wrong or missing centreline         Check that both points are in the
                                      same connected airway mask, confirm
                                      route direction, then review
                                      SlicerVMTK/Python errors.

  `no node was found`                 Check the exact Slicer node names
                                      and capitalisation expected by the
                                      script.

  Camera moves backwards              Confirm route point order and the
                                      start/end node assignments.

  Camera rolls in the wrong direction Check the camera roll convention
                                      and use the final verified script
                                      rather than an older version.

  CT and virtual endoscopic view      Stop at the last unambiguous
  disagree                            branch. Return to the CT and
                                      document the mismatch.

  Video export fails                  Save the scene, check the local
                                      video/FFmpeg setup if required by
                                      the script, then record again.

  Python error appears                Read the full error and determine
                                      whether the requested output was
                                      created correctly. Do not simply
                                      ignore the error. If using AI for
                                      debugging, paste the error text
                                      without patient data.
  -----------------------------------------------------------------------

------------------------------------------------------------------------

## 14. Final quality-control checklist

Before calling a case complete, confirm all of the following:

-   [ ] The correct CT series was selected.
-   [ ] The lesion marker is correct in axial, coronal and sagittal
    views.
-   [ ] The segmentation is derived from the chosen CT.
-   [ ] Every airway on the proposed route can be traced back to the CT
    or to a documented clinician correction.
-   [ ] No unsupported distal airway has been invented.
-   [ ] `RouteStart` is in the trachea.
-   [ ] `RouteTarget` is inside the final CT-supported segmented airway.
-   [ ] Start and target are connected within the verified airway
    segmentation.
-   [ ] The centreline stays within the intended airway route.
-   [ ] Each bifurcation has been checked against CT.
-   [ ] The virtual bronchoscopy has been reviewed from start to finish.
-   [ ] Any mismatch or uncertainty has been recorded.
-   [ ] Manual edits have been recorded.
-   [ ] Software, extension and script versions have been recorded.
-   [ ] The output is labelled as pre-procedural CT-based planning, not
    live navigation.

------------------------------------------------------------------------

## 15. What is automated and what still requires a person

The current workflow can automate much of the repetitive processing, but
it is not fully automatic.

**Automated or script-assisted:** - CT series ranking; - airway
segmentation; - 3D display; - centreline creation; - camera movement
along the route; - virtual bronchoscopy generation.

**Clinician input still required:** - choosing and checking the correct
CT; - marking the lesion; - checking the airway segmentation; -
selecting the final airway near the lesion; - placing `RouteStart` and
`RouteTarget`; - correcting genuine segmentation errors; - checking each
branch decision; - deciding when CT anatomy is too uncertain to
continue.

The aim is not to replace bronchial branch reading. It is to convert
CT-supported branch reading into a reproducible 3D and virtual
endoscopic roadmap.

------------------------------------------------------------------------

## 16. Reference documentation and background

Official software documentation:

-   3D Slicer DICOM guide:
    https://slicer.readthedocs.io/en/latest/user_guide/modules/dicom.html
-   3D Slicer Markups guide:
    https://slicer.readthedocs.io/en/latest/user_guide/modules/markups.html
-   3D Slicer Endoscopy guide:
    https://slicer.readthedocs.io/en/latest/user_guide/modules/endoscopy.html
-   SlicerTotalSegmentator extension:
    https://github.com/lassoan/SlicerTotalSegmentator
-   TotalSegmentator: https://github.com/wasserth/TotalSegmentator
-   SlicerVMTK: https://github.com/vmtk/SlicerExtension-VMTK

Clinical background for manual airway mapping and bronchial branch
tracing:

-   Kho SS, Nyanti LE, Chai CS, et al. Feasibility of manual bronchial
    branch reading technique in navigating conventional rEBUS
    bronchoscopy in the evaluation of peripheral pulmonary lesion. *Clin
    Respir J*. 2021;15:595-603.
-   Quah P, Lee P. Peripheral bronchoscopy: back to basics. *Curr Opin
    Pulm Med*. 2025. doi:10.1097/MCP.0000000000001229.
-   Kurimoto N, Morita K. *Bronchial Branch Tracing*. Singapore:
    Springer; 2020.

Software interfaces and model behaviour can change. For reproducibility,
record the exact versions used for each case and retain the final
scripts with the project.
