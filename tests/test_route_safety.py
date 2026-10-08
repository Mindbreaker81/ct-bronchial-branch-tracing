"""
Tests for the pure helper functions of the BBT Slicer scripts.

The scripts are standalone files meant to be executed inside the 3D Slicer
Python console, so they cannot be imported outside Slicer. These tests
extract only the named top-level functions from the script source and run
them with numpy, without Slicer.

Run with:  python -m pytest tests
"""

import ast
import math
import re
from pathlib import Path

import numpy as np
import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"

ROUTE_SCRIPT = SCRIPTS / "03_BBT_Continuous_Smooth_Route_v4_0.py"
BRONCH_SCRIPT = SCRIPTS / "04_BBT_Virtual_Bronch_v7_0.py"


def load_functions(path, names):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    selected = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
    missing = set(names) - {node.name for node in selected}
    assert not missing, f"functions not found in {path.name}: {missing}"
    namespace = {"np": np, "re": re, "math": math}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), "exec"), namespace)
    return namespace


route = load_functions(
    ROUTE_SCRIPT,
    ["minimumDistanceMM", "makesProgressTowardsTarget", "pointsOutsideMask"],
)
bronch = load_functions(BRONCH_SCRIPT, ["parseClockPosition"])


@pytest.mark.parametrize("path", sorted(SCRIPTS.glob("*.py")))
def test_scripts_compile(path):
    compile(path.read_text(encoding="utf-8"), str(path), "exec")


def test_minimum_distance_uses_physical_spacing():
    coords = np.array([[0, 0, 0], [2, 0, 0]])
    # Point is 1 voxel away along k (z) from [2,0,0]; z spacing is 2.5 mm.
    distance = route["minimumDistanceMM"](coords, [3, 0, 0], [2.5, 0.7, 0.7])
    assert distance == pytest.approx(2.5)


# Start-connected airway: a straight segment along k from 0 to 10.
START = np.array([[k, 20, 20] for k in range(0, 11)])
TARGET = np.array([30, 20, 20])
SPACING = np.array([1.0, 1.0, 1.0])
START_TO_TARGET = route["minimumDistanceMM"](START, TARGET, SPACING)


def test_start_front_distance():
    assert START_TO_TARGET == pytest.approx(20.0)


def test_fragment_ahead_of_route_front_is_bridgeable():
    fragment = np.array([[k, 20, 20] for k in range(14, 18)])
    assert route["makesProgressTowardsTarget"](fragment, START_TO_TARGET, TARGET, SPACING)


def test_target_component_is_always_bridgeable():
    target_component = np.array([[k, 20, 20] for k in range(25, 35)])
    assert route["makesProgressTowardsTarget"](target_component, START_TO_TARGET, TARGET, SPACING)


def test_speck_behind_route_front_is_rejected():
    speck = np.array([[2, 23, 20]])
    assert not route["makesProgressTowardsTarget"](speck, START_TO_TARGET, TARGET, SPACING)


def test_lateral_speck_beside_route_front_is_rejected():
    # Same depth as the front but 3 mm aside: farther from the target.
    speck = np.array([[10, 23, 20]])
    assert not route["makesProgressTowardsTarget"](speck, START_TO_TARGET, TARGET, SPACING)


def test_points_outside_mask_counts_out_of_bounds_and_background():
    mask = np.zeros((5, 5, 5), dtype=bool)
    mask[1:4, 1:4, 1:4] = True
    points = [
        [2.0, 2.0, 2.0],   # inside
        [2.4, 2.4, 2.4],   # rounds to inside
        [0.0, 0.0, 0.0],   # background
        [9.0, 2.0, 2.0],   # out of bounds
        [-1.0, 2.0, 2.0],  # negative index must not wrap around
    ]
    assert route["pointsOutsideMask"](points, mask) == 3


@pytest.mark.parametrize(
    "description, expected",
    [
        ("3", 3),
        ("12", 12),
        ("3 o'clock", 3),
        ("3, 3", 3),
        ("", None),
        (None, None),
        ("no value", None),
        ("13", None),
        ("0", None),
        ("branch 6, 3 o'clock", None),
    ],
)
def test_parse_clock_position(description, expected):
    assert bronch["parseClockPosition"](description) == expected
