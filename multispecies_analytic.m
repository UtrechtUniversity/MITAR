%% Multispecies gLV without and with plasmids %%

%% Preparations for model without plasmids
reset(symengine)
clear
clc

% Define variables and equations
syms a11 a12 a21 a22 R1 R2 r1 r2
assume([a11, a22] < 0)
intmat = [[a11, a12]; [a21, a22]];
mu = [r1; r2];
abun = [R1; R2];
fun2 = diag(abun) * (mu + intmat * abun) % The 2 denotes it is intended for use with 2 species


%% Calculate equilibrium abundances
Eq = solve(fun2 == 0, [R1, R2], 'ReturnConditions', true);
[Eq.R1, Eq.R2]

% Select feasible equilibrium, i.e., both R1 and R2 are larger than zero. 'find' returns indices of
% nonzero elements, using intersect to get indices of values that are simultaneously nonzero for R1
% and R2
indexfeasible = intersect(find(Eq.R1), find(Eq.R2))
Eq1 = Eq.R1(indexfeasible)
Eq2 = Eq.R2(indexfeasible)
assumeAlso(Eq.conditions(indexfeasible))

% Calculate growth rates required to obtain an equilibrium with the given abundances and
% interactions
sol = solve(fun2 == 0, [r1, r2]);
mu1 = sol.r1
mu2 = sol.r2


%% Stability analysis
% Commented this section out, since it is the same as the ecological
% stability in the model with plasmids analysed below

% % Calculate Jacobian matrix and corresponding eigenvalues at coexistence equilibrium
% JacgLV = simplify(jacobian([fun2], [R1, R2]), 'Steps', 200) % general jacobian
% % Jacobian at coexistence equilibrium
% JacgLVReduced = simplify(subs(JacgLV, [R1, R2], [Eq1, Eq2]), 'Steps', 200)
% % Eigenvalues of jacobian at coexistence equilibrium
% eigValsReduced = simplify(eig(JacgLVReduced), 'Steps', 200)

% % Jacobian at coexistence equilibrium considering used parameterization
% JacgLVReducedParm = simplify(subs(JacgLVReduced, [r1, r2], [mu1, mu2]), 'Steps', 200)
% eigValsJacgLVReducedParm = simplify(eig(JacgLVReducedParm), 'Steps', 200)
% % See below on ecological stability for derivation of stability criterion a11*a22 < a12*a21

%% Preparations model with plasmids
% Not clearing memory and symbolic toolbox because the equilibrium abundances calculated above
% will be used as abundances at the plasmid-free equilibrium

% Define additional variables and equations
syms P1 P2
syms c1 c2 g11 g12 g21 g22
assume([c1, c2, g11, g12, g21, g22] > 0)

dR1 = R1 * (r1        + a11*(R1 + P1) + a12*(R2 + P2)) - R1*(g11*P1 + g12*P2)
dR2 = R2 * (r2        + a21*(R1 + P1) + a22*(R2 + P2)) - R2*(g21*P1 + g22*P2)    
dP1 = P1 * (r1 - c1 + a11*(R1 + P1) + a12*(R2 + P2)) + R1*(g11*P1 + g12*P2)
dP2 = P2 * (r2 - c2 + a21*(R1 + P1) + a22*(R2 + P2)) + R2*(g21*P1 + g22*P2)


%% Analysis of general stability
% Calculate Jacobian matrix and corresponding eigenvalues at plasmid-free equilibrium
% general jacobian
JacGen = simplify(jacobian([dR1, dR2, dP1, dP2], [R1, R2, P1, P2]), 'Steps', 200)
% Jacobian if no plasmid-bearing bacteria are present
JacPfree = simplify(subs(JacGen, [P1, P2], [0, 0]), 'Steps', 200)
% Jacobian at plasmid-free equilibrium
JacPfreeEq = simplify(subs(JacPfree, [R1, R2], [Eq1, Eq2]), 'Steps', 200)
% Jacobian at plasmid-free equilibrium considering used parameterization
JacPfreeEqParm = simplify(subs(JacPfreeEq, [r1, r2], [mu1, mu2]), 'Steps', 200)

% Next lines commented out since it takes long
% % Eigenvalues of jacobian at plasmid-free equilibrium
% eigValsPfreeEq = simplify(eig(JacPfreeEq), 'Steps', 200)
% % Eigenvalues of jacobian at plasmid-free equilibrium considering used parameterization
% eigValsPfreeEqParm = simplify(eig(JacPfreeEqParm), 'Steps', 200)

%% Analysis of ecological stability
if(unique(JacPfree(3:4, 1:2)) ~= 0) 
    error("The lower part of JacPfree should be only zeros!")
end

JacPfreeEcol = JacPfree(1:2, 1:2) % At plasmid-free equilibrium
% At plasmid-free equilibrium considering used parameterization
JacPfreeEqParmEcol = JacPfreeEqParm(1:2, 1:2)
eigValsJacPfreeEqParmEcol = simplify(eig(JacPfreeEqParmEcol), 'Steps', 200)

% Rewrite eigenvalue to get a stability criterion. Equilibrium is ecologically stable if both
% eigenvalues are negative, so the maximum eigenvalue is the relevant one to derive
% the criterion. Starting with LHS being the eigenvalue and RHS = 0, so the equilibrium
% is ecologically stable if LHS < RHS (but inequality changes during the derivation, see
% below)
if max(eigValsJacPfreeEqParmEcol, [], 'all') == eigValsJacPfreeEqParmEcol(1)
    IndexLargestEigval = 1
    disp("FIRST eigenvalue is largest")
else
    IndexLargestEigval = 2
    disp("SECOND eigenvalue is largest")
end

releigValLHS = eigValsJacPfreeEqParmEcol(IndexLargestEigval)
releigValRHS = 0

% Subtract term with the root from both sides
[coefsLHS, termsLHS] = coeffs(releigValLHS);
releigValLHS = releigValLHS - coefsLHS(3) * termsLHS(3)
releigValRHS = releigValRHS - coefsLHS(3) * termsLHS(3)
% stable if LHS < RHS

% Multiply with -2 to get rid of the fractions and make both sides positive (reverse inequality
% because multiplying with a negative number), and square both sides.
releigValLHSsq = expand((-2 * releigValLHS)^2)
releigValRHSsq = expand((-2 * releigValRHS)^2)
% stable if LHS > RHS (reversed inequality!)

releigValLHSsq = releigValLHSsq - R1^2*a11^2
releigValRHSsq = releigValRHSsq - R1^2*a11^2
% stable if LHS > RHS

releigValLHSsq = releigValLHSsq - R2^2*a22^2
releigValRHSsq = releigValRHSsq - R2^2*a22^2
% stable if LHS > RHS

releigValLHSsq = releigValLHSsq + 2*R1*R2*a11*a22
releigValRHSsq = releigValRHSsq + 2*R1*R2*a11*a22
% stable if LHS > RHS

releigValLHSsq = releigValLHSsq / (4*R1*R2) % a11*a22
releigValRHSsq = releigValRHSsq / (4*R1*R2) % a12*a21
% stable if LHS > RHS: a11*a22 > a12*a21; unstable if  LHS < RHS: a11*a22 < a12*a21

% Determining stability criterion for ecological stability with the determinant-trace method
simplify(det(JacPfreeEqParmEcol), 'Steps', 100) % product of the eigenvalues, > 0 for stable point
% R1*R2*(a11*a22 - a12*a21) > 0 -> a11*a22 > a12*a21 
trace(JacPfreeEqParmEcol) % sum of the eigenvalues, < 0 for stable point
% R1*a11 + R2*a22 -> always < 0 because (a11, a22) < 0 and (R1, R2) > 0
% So same result as from the derivation above:
% stable if a11*a22 > a12*a21; unstable if LHS < RHS: a11*a22 < a12*a21

%% Analysis of epidemiological stability: full model
JacPfreeEpi = JacPfree(3:4, 3:4) % At plasmid-free equilibrium
% At plasmid-free equilibrium considering used parameterization
JacPfreeEqParmEpi = JacPfreeEqParm(3:4, 3:4)
eigValsPfreeEqParmEpi = simplify(eig(JacPfreeEqParmEpi), 'Steps', 200)

% Determining stability criterion for epidemiological stability with the determinant-trace method
traceEpi = trace(JacPfreeEqParmEpi)
% R1*g11 - c2 - c1 + R2*g22
detEpi = det(JacPfreeEqParmEpi)
% c1*c2 - R1*c2*g11 - R2*c1*g22 + R1*R2*g11*g22 - R1*R2*g12*g21
discrEpi = traceEpi^2 - 4 * detEpi
% (c1 + c2 - R1*g11 - R2*g22)^2 - 4*c1*c2 + 4*R1*c2*g11 + 4*R2*c1*g22 - 4*R1*R2*g11*g22 +
% 4*R1*R2*g12*g21
eigEpi1 = (traceEpi + sqrt(discrEpi)) / 2
eigEpi2 = (traceEpi - sqrt(discrEpi)) / 2
% Stable if det > 0 and tr < 0

simplify(eigValsPfreeEqParmEpi(1) - eigEpi2, 'Steps', 200) % 0
simplify(eigValsPfreeEqParmEpi(2) - eigEpi1, 'Steps', 200) % 0

% (R1*g11)/2 - c2/2 - c1/2 + (R2*g22)/2 +/- (R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 -
% 2*R1*c1*g11 + 2*R1*c2*g11 + R2^2*g22^2 + 2*R2*c1*g22 - 2*R2*c2*g22 + c1^2 - 2*c1*c2 +
% c2^2)^(1/2)/2 < 0

% Both eigenvalues have to be negative to ensure stability, so the relevant eigenvalue to determine
% stability is the one where the square root is added to the other part
% (R1*g11)/2 - c2/2 - c1/2 + (R2*g22)/2 + (R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 - 
% 2*R1*c1*g11 + 2*R1*c2*g11 + R2^2*g22^2 + 2*R2*c1*g22 - 2*R2*c2*g22 + c1^2 - 2*c1*c2 + 
% c2^2)^(1/2)/2 < 0

% (R1*g11)/2 - c2/2 - c1/2 + (R2*g22)/2 < -1*(R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 -
% 2*R1*c1*g11 + 2*R1*c2*g11 + R2^2*g22^2 + 2*R2*c1*g22 - 2*R2*c2*g22 + c1^2 - 2*c1*c2 + 
% c2^2)^(1/2)/2

% R1*g11 - c2 - c1 + R2*g22 < -1*(R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 -
% 2*R1*c1*g11 + 2*R1*c2*g11 + R2^2*g22^2 + 2*R2*c1*g22 - 2*R2*c2*g22 + c1^2 - 2*c1*c2 +
% c2^2)^(1/2)
R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 - 2*R1*c1*g11 + 2*R1*c2*g11 + R2^2*g22^2 + 2*R2*c1*g22 - 2*R2*c2*g22 + c1^2 - 2*c1*c2 + c2^2 % This line should be removed?
% R1*g11 - c2 - c1 + R2*g22 < -1*(((R1*g11 - c1) - (R2*g22 - c2))^2 + 4*g12*g21*R1*R2)^(1/2)

% This last rearrangement shows the part inside the root is always positive, so the criterion can be
% simplified by squaring both sides without the need to reverse the sign of the inequality
releigValepiLHSsq = (R1*g11 - c2 - c1 + R2*g22)^2
releigValepiRHSsq = ((R1*g11 - c1) - (R2*g22 - c2))^2 + 4*g12*g21*R1*R2
% stable if LHSsq < RHSsq, so stable if LHSsq - RHSsq < 0

LHSsq = simplify(expand(releigValepiLHSsq) - expand(releigValepiRHSsq), 'Steps', 200)
% 4*c1*c2 - 4*R1*c2*g11 - 4*R2*c1*g22 + 4*R1*R2*g11*g22 - 4*R1*R2*g12*g21 < 0
RHSsq = 0

LHSsq = LHSsq / 4
RHSsq = RHSsq / 4
% c1*c2 - R1*c2*g11 - R2*c1*g22 + R1*R2*g11*g22 - R1*R2*g12*g21
% stable if LHSsq < RHSsq

LHSsq = LHSsq - c1*c2
RHSsq = RHSsq - c1*c2
% stable if LHSsq < RHSsq

LHSsq = -1*LHSsq
RHSsq = -1*RHSsq
% stable if LHSsq > RHSsq (reversed inequality!)
% stable if R1*c2*g11 + R2*c1*g22 - R1*R2*g11*g22 + R1*R2*g12*g21 > c1*c2

%% Analysis of epidemiological stability: c1 = c2 = c
syms c
dP1c = simplify(subs(dP1, [c1, c2], [c, c]),  'Steps', 200)
dP2c = simplify(subs(dP2, [c1, c2], [c, c]),  'Steps', 200)

JacGenc = simplify(jacobian([dR1, dR2, dP1c, dP2c], [R1, R2, P1, P2]), 'Steps', 200) % general jacobian
% Jacobian if no plasmid-bearing bacteria are present
JacPfreec = simplify(subs(JacGenc, [P1, P2], [0, 0]), 'Steps', 200)
% Jacobian at plasmid-free equilibrium
JacPfreeEqc = simplify(subs(JacPfreec, [R1, R2], [Eq1, Eq2]), 'Steps', 200)
% Jacobian at plasmid-free equilibrium considering used parameterization
JacPfreeEqParmc = simplify(subs(JacPfreeEqc, [r1, r2], [mu1, mu2]), 'Steps', 200)
% At plasmid-free equilibrium considering used parameterization
JacPfreeEqParmEpic = JacPfreeEqParmc(3:4, 3:4)
eigValsPfreeEqParmEpic = simplify(eig(JacPfreeEqParmEpic), 'Steps', 200)
% (R1*g11)/2 - (R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 + R2^2*g22^2)^(1/2)/2 - c + (R2*g22)/2
% (R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 + R2^2*g22^2)^(1/2)/2 - c + (R1*g11)/2 + (R2*g22)/2

% Reordering terms and using +/- notation
% (R1*g11)/2 + (R2*g22)/2 - c +/- (R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 + R2^2*g22^2)^(1/2)/2 

% Rewriting into stability crterion in terms of costs
% (R1*g11)/2 + (R2*g22)/2 - c +/- (R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 + R2^2*g22^2)^(1/2)/2 < 0
% (R1*g11)/2 + (R2*g22)/2 - c + (R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 + R2^2*g22^2)^(1/2)/2 < 0
% (R1*g11)/2 + (R2*g22)/2 + (R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 + R2^2*g22^2)^(1/2)/2 < c
% (R1*g11) + (R2*g22) + (R1^2*g11^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2 + R2^2*g22^2)^(1/2) < 2c
% (R1*g11) + (R2*g22) + (R1^2*g11^2 + R2^2*g22^2 - 2*R1*R2*g11*g22 + 4*g12*g21*R1*R2)^(1/2) < 2c
% (R1*g11) + (R2*g22) + ((R1*g11 - R2*g22)^2 + 4*g12*g21*R1*R2)^(1/2) < 2c
% 2c > R1*g11 + R2*g22 + sqrt((R1*g11 - R2*g22)^2 + 4*g12*g21*R1*R2)

% Rewriting into stability crterion in terms of conjugation rates
% 2c > R1*g11 + R2*g22 + sqrt((R1*g11 - R2*g22)^2 + 4*g12*g21*R1*R2)
% - R1*g11 - R2*g22 > sqrt((R1*g11 - R2*g22)^2 + 4*g12*g21*R1*R2) - 2c
% R1*g11 + R2*g22 < 2c - sqrt((R1*g11 - R2*g22)^2 + 4*g12*g21*R1*R2)


%% Analysis of epidemiological stability if g12 = 0, g21 = 0
dP1galt = simplify(subs(dP1, [g12, g21], [0, 0]),  'Steps', 200)
dP2galt = simplify(subs(dP2, [g12, g21], [0, 0]),  'Steps', 200)

JacGengalt = simplify(jacobian([dR1, dR2, dP1galt, dP2galt], [R1, R2, P1, P2]), 'Steps', 200) % general jacobian
% Jacobian if no plasmid-bearing bacteria are present
JacPfreegalt = simplify(subs(JacGengalt, [P1, P2], [0, 0]), 'Steps', 200)
% Jacobian at plasmid-free equilibrium
JacPfreeEqgalt = simplify(subs(JacPfreegalt, [R1, R2], [Eq1, Eq2]), 'Steps', 200)
% Jacobian at plasmid-free equilibrium considering used parameterization
JacPfreeEqParmgalt = simplify(subs(JacPfreeEqgalt, [r1, r2], [mu1, mu2]), 'Steps', 200)
% At plasmid-free equilibrium considering used parameterization
JacPfreeEqParmEpigalt = JacPfreeEqParmgalt(3:4, 3:4)
eigValsPfreeEqParmEpigalt = simplify(eig(JacPfreeEqParmEpigalt), 'Steps', 200)
% R1*g11 - c1
% R2*g22 - c2

% Rewriting into stability crterion in terms of costs
% R1*g11 - c1 < 0 and R2*g22 - c2 < 0
% R1*g11 < c1 and R2*g22 < c2
% c1 > R1*g11 and c2 > R2*g22

% Rewriting into stability crterion in terms of conjugation rates
% c1 > R1*g11 and c2 > R2*g22
% R1*g11 < c1 and R2*g22 < c2


%% Analysis of epidemiological stability if c1 = c2 = c and g12, g21 = 0
dP1cgalt = simplify(subs(dP1, [c1, c2, g12, g21], [c, c, 0, 0]),  'Steps', 200)
dP2cgalt = simplify(subs(dP2, [c1, c2, g12, g21], [c, c, 0, 0]),  'Steps', 200)

% general jacobian
JacGencgalt = simplify(jacobian([dR1, dR2, dP1cgalt, dP2cgalt], [R1, R2, P1, P2]), 'Steps', 200)
% Jacobian if no plasmid-bearing bacteria are present
JacPfreecgalt = simplify(subs(JacGencgalt, [P1, P2], [0, 0]), 'Steps', 200)
% Jacobian at plasmid-free equilibrium
JacPfreeEqcgalt = simplify(subs(JacPfreecgalt, [R1, R2], [Eq1, Eq2]), 'Steps', 200)
% Jacobian at plasmid-free equilibrium considering used parameterization
JacPfreeEqParmcgalt = simplify(subs(JacPfreeEqcgalt, [r1, r2], [mu1, mu2]), 'Steps', 200)
% At plasmid-free equilibrium considering used parameterization
JacPfreeEqParmEpicgalt = JacPfreeEqParmcgalt(3:4, 3:4)
eigValsPfreeEqParmEpicgalt = simplify(eig(JacPfreeEqParmEpicgalt), 'Steps', 200)
% R1*g11 - c
% R2*g22 - c

% Rewriting into stability crterion in terms of costs
% R1*g11 - c < 0 and R2*g22 - c < 0
% R1*g11 < c and R2*g22 < c
% max(R1*g11, R2*g22) < c
% c > max(R1*g11, R2*g22)

% Rewriting into stability crterion  in terms of conjugation rates
% c > max(R1*g11, R2*g22)
% max(R1*g11, R2*g22) < c


%% Analysis of submatrix D
JacPfreeD = JacPfree(1:2,  3:4) % At plasmid-free equilibrium
% At plasmid-free equilibrium considering used parameterization
JacPfreeEqParmD = JacPfreeEqParm(1:2,  3:4)
eigValsJacPfreeEqParmD = simplify(eig(JacPfreeEqParmD), 'Steps', 200)


%% Find the intersection of the bifurcation lines for the two-species and sixteen-species model
% To find the intraspecies conjugation rate for which the bifurcation liness cross.
syms g
assumeAlso(g > 0)

inv_crit = R1*g11 + R2*g22 + sqrt((R1*g11 - R2*g22)^2 + 4*g12*g21*R1*R2)
inv_crit_2sp = simplify(subs(inv_crit, [g22, g12, g21], [g, g, g]), 'Steps', 200)
% The factors 7.5 / 8 and 2.5 / 2 are based on the abundances of the initially plasmid-bearing
% species and the initially plasmid-free species in the two-species and sixteen-species models.
inv_crit_16sp = simplify(subs(inv_crit, [R1, R2, g22, g12, g21], [R1 * 7.5 / 8, R2 * 2.5 / 2, g, g, g]), 'Steps', 200)

sol = solve(inv_crit_2sp == inv_crit_16sp, [g11], 'ReturnConditions', true);
intersect_g11 = simplify(sol.g11, 'Steps', 200)
sol.conditions

simplify(subs(intersect_g11, [R1, R2], [8 * 10^10, 2 * 10^10]), 'Steps', 200) % g
% So the bifurcation lines cross where g11 == g, as is indeed shown in the bifurcation-like plot
% based on numerical simulations.

% Now solve the invasion criterion for case where all conjugation rates are equal to each other.
inv_crit_all_g = simplify(subs(inv_crit, [g11, g22, g12, g21], [g, g, g, g]), 'Steps', 200)
% g*(R1 + R2 + ((R1 + R2)^2)^(1/2))


%% Solve for g11
asol1 = solve(eigEpi1 == 0, [g11], 'ReturnConditions', true)
simplify(asol1.g11, 'Steps', 100) % c1/R1 - (R2*g12*g21)/(c2 - R2*g22)

asol2 = solve(eigEpi2 == 0, [g11], 'ReturnConditions', true)
simplify(asol2.g11, 'Steps', 100) % c1/R1 - (R2*g12*g21)/(c2 - R2*g22)
