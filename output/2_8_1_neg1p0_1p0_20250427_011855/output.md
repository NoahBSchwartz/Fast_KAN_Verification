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
Defined 565 base constraints.

=== MIP Problem Complexity ===
Network Dimensions: [2, 8, 1]
Total Variables: 168
  - Continuous Variables: 35
  - Binary Variables: 133
Total Constraints: 565
Total Segments: 133
  - Layer 0: 71 segments (avg: 4.44, max: 8)
  - Layer 1: 62 segments (avg: 7.75, max: 9)
Indicator Variables / Total Variables: 79.17%
============================


--- Optimizing for Output Dimension 0 ---
Solving for Min Output 0 using Gurobi...
Min Output 0 found: -0.645010
Solving for Max Output 0 using Gurobi...
Max Output 0 found: 0.064409
Lipschitz Verification Time: 0.1757 seconds
Lipschitz Total Time (conversion + verification): 0.1780 seconds

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
Min Output 0 found: -0.566771
Solving for Max Output 0 using Gurobi...
Max Output 0 found: -0.024932
Vanilla Verification Time: 0.2796 seconds
Vanilla Total Time (conversion + verification): 0.2813 seconds

--------------------------------------------------
VALIDATING MIP BOUNDS
--------------------------------------------------

Generating 5000 validation samples...
Running KAN model on validation samples...

Validating vanilla MIP bounds...
Output 0:
  MIP bounds:    [-0.566771, -0.024932]
  Sampled bounds: [-0.419116, -0.193329]
  ✓ Bounds valid

✓ All vanilla MIP bounds are valid

Validating Lipschitz MIP bounds...
Output 0:
  MIP bounds:    [-0.645010, 0.064409]
  Sampled bounds: [-0.419116, -0.193329]
  ✓ Bounds valid

✓ All Lipschitz MIP bounds are valid

VERIFICATION PASSED FOR BOTH METHODS
