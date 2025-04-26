================================================================================
VERIFICATION SETTINGS
================================================================================
Input Bounds: [(-0.3, 0.3), (-0.3, 0.3)]
Time Limit (seconds): 6000
Validate Bounds: True
Number of Validation Samples: 5000
Absolute Tolerance: 0.001
Relative Tolerance: 0.01

--------------------------------------------------
VERIFICATION COMPARISON
--------------------------------------------------

Converting vanilla segments to verification format...
Processing Layer 0 segments (16 entries)...
Processing Layer 1 segments (8 entries)...
Segment conversion complete.

Running MIP verification for vanilla segments...
Setting up MIP for KAN [2, 8, 1]...
Input bounds: [(-0.3, 0.3), (-0.3, 0.3)]
Defined 321 base constraints.

=== MIP Problem Complexity ===
Network Dimensions: [2, 8, 1]
Total Variables: 107
  - Continuous Variables: 35
  - Binary Variables: 72
Total Constraints: 321
Total Segments: 72
  - Layer 0: 48 segments (avg: 3.00, max: 3)
  - Layer 1: 24 segments (avg: 3.00, max: 3)
Indicator Variables / Total Variables: 67.29%
============================


--- Optimizing for Output Dimension 0 ---
Solving for Min Output 0 using Gurobi...
Min Output 0 found: -1.177760
Solving for Max Output 0 using Gurobi...
Max Output 0 found: 0.328119
Vanilla Verification Time: 0.0781 seconds
Vanilla Total Time (conversion + verification): 0.0808 seconds

--------------------------------------------------
VALIDATING MIP BOUNDS
--------------------------------------------------

Generating 5000 validation samples...
Running KAN model on validation samples...

Validating vanilla MIP bounds...
Output 0:
  MIP bounds:    [-1.177760, 0.328119]
  Sampled bounds: [-0.656549, -0.582623]
  ✓ Bounds valid

✓ All vanilla MIP bounds are valid

VERIFICATION PASSED FOR BOTH METHODS
