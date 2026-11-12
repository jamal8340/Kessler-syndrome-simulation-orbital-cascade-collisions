# Review of the Supplied Unity Project

## Scope

The following files were inspected from `KesslerSimulation-main.zip`: -
`Assets/SatelliteSpawner.cs`; - `Assets/OrbitCameraController.cs`; -
`Assets/Scenes/SampleScene.unity`; - `Assets/Satelite.prefab`; -
`Assets/StreamingAssets/satellites_state.json`; -
`Packages/manifest.json`; - `ProjectSettings/ProjectVersion.txt`.

The project uses Unity `6000.6.7f1` and URP. The scene contains
`Main Camera`, `Directional Light`, `SatelliteSpawner`, `Global Volume`,
`Manager` and `Earth`.

## What works and should be retained

### Visualization and controls

-   3D scene with Earth and its texture.
-   Satellite prefab.
-   Orbit camera with rotation and zoom.
-   Simulation speed controls and pause.
-   JSON array loading from `StreamingAssets`.
-   Object instantiation and category-based coloring.
-   On-screen statistics panel.
-   Disabled prefab colliders, avoiding Unity PhysX processing thousands
    of orbital objects.

### Orbital motion

`SatelliteSpawner.Update()` stores positions and velocities in km and
km/s. For each object it computes central Earth gravity:

`a = -mu * r / |r|^3`

It then updates velocity and position:

`v[n+1] = v[n] + a[n] * dt`

`r[n+1] = r[n] + v[n+1] * dt`

This is semi-implicit (symplectic) Euler integration, not SGP4
propagation. It is simpler and often has better energy behaviour than
explicit Euler, but it is still a first-order method. Its error depends
on the integration step.

## Main physical and numerical issues

### 1. The integration step can be too large

The code caps frame delta time and then divides it into at most 10
substeps. At high time acceleration, each substep can still be several
seconds. For fast orbital objects this can create significant position
error and missed close approaches.

### 2. The spatial hash does not inspect neighbouring cells

`CheckCollisions()` compares an object only with earlier objects
assigned to exactly the same spatial hash key. It does not inspect the
26 neighbouring cells in 3D.

Two objects closer than the collision threshold can lie on opposite
sides of a cell boundary and never be compared. This creates important
false negatives.

### 3. Detection uses instantaneous separation, not TCA

The code checks distance only at the instant `CheckCollisions()` runs.
It does not analyse relative motion between previous and current steps
and does not calculate Time of Closest Approach (TCA).

A fast encounter can occur and end between checks.

### 4. A 2 km threshold is not a physical collision radius

`collisionDistanceKm = 2.0` treats objects within 2 km as colliding.
This can be used as a broad conjunction-screening threshold, but not as
a physical contact condition. Physical contact requires object
dimensions and a separation comparable to the sum of their radii.

### 5. Minimum relative velocity filter

`minRelativeVelocityKmS = 1.0` rejects pairs with lower relative speed.
This may suppress some false alarms, but it is not a general physical
condition for collision. A slower encounter can also result in contact.

### 6. Collision handling does not model fragmentation

`HandleCollision()` marks two existing objects as debris.
`DestroySatelliteIntoDebris()` adds a random velocity impulse to an
existing object.

The code does not: - create new fragments; - distribute mass; -
calculate collision energy; - conserve system momentum; - use an
empirical fragmentation model; - propagate fragment size and velocity
distributions.

This is a visual event placeholder, not a physical cascade model.

### 7. The C-key test is artificial

`TriggerTestCollision()` selects a random active object in a specified
LEO altitude range and changes its status. It does not detect an actual
encounter between two objects. It is useful for UI demonstration but
must not be treated as a scientific result.

### 8. Reference frame and axis mapping need verification

SGP4 returns state vectors in TEME. Unity copies `x, y, z` directly
without an explicit transformation. Unity commonly uses Y as its up
axis, while geocentric orbital data commonly use Z as the polar axis in
the chosen reference frame.

The Earth model's orientation must be checked and an explicit axis
transformation defined. Direct XYZ assignment should not be assumed to
be geographically correct.

Rotating Earth around `Vector3.up` does not by itself implement a
correct TEME-to-ECEF transformation.

### 9. Input state and epoch

Unity loads a static `satellites_state.json`, but the data structure
used by `SatelliteSpawner` has no epoch timestamp or reference-frame
metadata. This makes reproducibility and comparison with Python
difficult.

### 10. The bundled Unity data differs from the latest Python output

`Assets/StreamingAssets/satellites_state.json` contains 18,579
objects: - 15,897 `active`; - 2,677 `debris`; - 5 `rocket_body`.

The latest reported Python run processed 18,479 objects. These are
different snapshots or processing states. One canonical data file should
be generated and used by both systems.

### 11. Numeric precision

Unity uses `Vector3` and `float` for positions in kilometres. This may
be adequate for LEO visualization, but precise close-approach
calculations should use `double`; convert to float only after
subtracting vectors and preparing render coordinates.

## Overall assessment

  -----------------------------------------------------------------------
  Area                                Assessment
  ----------------------------------- -----------------------------------
  3D scene and visualization          Good prototype foundation

  JSON catalog loading                Works, but lacks epoch/frame
                                      metadata

  Central gravity motion              Simplified two-body model

  Numerical integration               Semi-implicit Euler; step-size
                                      tests required

  Collision screening                 Incomplete; significant false
                                      negatives possible

  Physical contact                    Not implemented

  Fragmentation                       Placeholder; no physical accounting

  Kessler cascade                     Event demonstration, not a
                                      validated model

  Python integration                  Not implemented yet
  -----------------------------------------------------------------------

## Recommended direction

Do not rewrite all of Unity. Retain: - scene, Earth, textures and
prefab; - camera, time controls and GUI; - object rendering; - visual
event indicators.

Gradually replace or subordinate the computation layer to Python backend
outputs: 1. Shared data format with NORAD ID, UTC epoch, reference frame
and units. 2. One canonical object snapshot. 3. Python
close-approach/TCA module. 4. Events and trajectories exported to Unity.
5. Fragmentation and cascade model only after that.

Unity may retain a simplified demonstration mode, but it should be
clearly labelled as visualization mode rather than scientific mode.
