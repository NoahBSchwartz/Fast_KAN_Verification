================================================================================
VERIFICATION SETTINGS
================================================================================
Input Bounds: [(-1.0, 1.0), (-1.0, 1.0)]
Time Limit (seconds): 6000
Validate Bounds: True
Number of Validation Samples: 5000
Absolute Tolerance: 0.001
Relative Tolerance: 0.01

--------------------------------------------------
VERIFICATION COMPARISON
--------------------------------------------------

Converting Lipschitz segments to verification format...
Processing Layer 0 segments (16 entries)...
Processing Layer 1 segments (8 entries)...
Segment conversion complete.

Running MIP verification for Lipschitz segments...
Setting up MIP for KAN [2, 8, 1]...
Input bounds: [(-1.0, 1.0), (-1.0, 1.0)]
Defined 557 base constraints.

=== MIP Problem Complexity ===
Network Dimensions: [2, 8, 1]
Total Variables: 166
  - Continuous Variables: 35
  - Binary Variables: 131
Total Constraints: 557
Total Segments: 131
  - Layer 0: 65 segments (avg: 4.06, max: 6)
  - Layer 1: 66 segments (avg: 8.25, max: 10)
Indicator Variables / Total Variables: 78.92%
============================


--- Optimizing for Output Dimension 0 ---
Solving for Min Output 0 using Gurobi...
Min Output 0 found: -0.377545
Solving for Max Output 0 using Gurobi...
Max Output 0 found: 0.397822
Lipschitz Verification Time: 0.1224 seconds
Lipschitz Total Time (conversion + verification): 0.1247 seconds

Converting vanilla segments to verification format...
Processing Layer 0 segments (16 entries)...
Processing Layer 1 segments (8 entries)...
Segment conversion complete.

Running MIP verification for vanilla segments...
Setting up MIP for KAN [2, 8, 1]...
Input bounds: [(-1.0, 1.0), (-1.0, 1.0)]
Defined 993 base constraints.

=== MIP Problem Complexity ===
Network Dimensions: [2, 8, 1]
Total Variables: 275
  - Continuous Variables: 35
  - Binary Variables: 240
Total Constraints: 993
Total Segments: 240
  - Layer 0: 160 segments (avg: 10.00, max: 10)
  - Layer 1: 80 segments (avg: 10.00, max: 10)
Indicator Variables / Total Variables: 87.27%
============================


--- Optimizing for Output Dimension 0 ---
Solving for Min Output 0 using Gurobi...
Min Output 0 found: -0.245412
Solving for Max Output 0 using Gurobi...
Max Output 0 found: 0.261836
Vanilla Verification Time: 0.3634 seconds
Vanilla Total Time (conversion + verification): 0.3649 seconds

--------------------------------------------------
VALIDATING MIP BOUNDS
--------------------------------------------------

Generating 5000 validation samples...
Running KAN model on validation samples...

Validating vanilla MIP bounds...
Output 0:
  MIP bounds:    [-0.245412, 0.261836]
  Sampled bounds: [-0.098442, 0.085463]
  ✓ Bounds valid

✓ All vanilla MIP bounds are valid

Validating Lipschitz MIP bounds...
Output 0:
  MIP bounds:    [-0.377545, 0.397822]
  Sampled bounds: [-0.098442, 0.085463]
  ✓ Bounds valid

✓ All Lipschitz MIP bounds are valid

VERIFICATION PASSED FOR BOTH METHODS
