# Changelog

## [Unreleased]

Derived from v1.0.0 (`3567e86`). Not yet reproduced inside 3D Slicer.

### Changed
- Route script 03: stop when `RouteStart` or `RouteTarget` is more than
  5 mm from the airway segment instead of snapping silently.
- Route script 03: bridge only to fragments closer to the target than the
  start-connected airway; isolated fragments behind or beside the route
  front can no longer consume bridges.
- Route script 03: centre slice views on each gap before the bridge
  authorisation dialog.
- Route script 03: constrain the smoothed route to the Start-to-Target
  connected component instead of the whole navigation mask.
- Virtual bronchoscopy script 04: use `BBT.RouteDirection` when present
  instead of the head-to-feet heuristic.
- Virtual bronchoscopy script 04: ignore `B<n>` markers more than 10 mm
  from the route, or whose description has no single clock value 1-12.

### Added
- Route script 03: count of final route points outside the route
  component, printed and stored in `BBT.PointsOutsideRouteComponent`.
- Route script 03: `BBT.RouteDirection` attribute on `Route_Centerline`.
- `tests/` with Slicer-free tests of the pure helper functions.
- Guide: built-in safety checks, route direction and bifurcation marker
  conventions.

## [1.0.0] - 2026-10-06

### Added
- Public repository structure
- Reproduction documentation
- Citation framework
- Contribution guidance
- Safety disclaimer
- Validation templates

### Required before public release
- Insert exact tested Python scripts
- Confirm Slicer/extension version manifest
- Add full reproduction checklist PDF
- Replace licence placeholder with complete official Apache-2.0 text
- Archive v1.0.0 with Zenodo and add DOI
