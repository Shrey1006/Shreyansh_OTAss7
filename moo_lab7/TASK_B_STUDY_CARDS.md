# Task B: Algorithm Study Cards

Use each card as a short viva prompt: explain the mechanism, then name the trade-off or implementation detail. The descriptions match the concepts represented in this project; compact variants are identified in the project README and results report.

| ID | Algorithm / family | Recall card | Trade-off / viva check |
|---:|---|---|---|
| 1 | Weighted Sum - classical | Minimise a weighted sum of normalised objectives; sweep weights to sample preferences. | A linear scalarisation can miss non-convex portions of a Pareto front. |
| 2 | Epsilon-Constraint - classical | Optimise one objective while treating other objectives as swept upper bounds. | Explain how epsilon values move the feasible region and why infeasible scalar runs need penalties. |
| 3 | Benson's Method - classical | Improve a feasible reference point in objective space while seeking feasible objective-wise improvement. | Results depend on the reference point; this implementation uses repeated scalar subproblems. |
| 4 | Tchebycheff Weight Metric - classical | Minimise the maximum weighted deviation from the ideal point. | Unlike a weighted sum, it can recover non-convex front regions; weight and scaling choices still matter. |
| 5 | Value Function - classical | Maximise a weighted Cobb-Douglas-style utility over normalised objective values. | Weights encode preferences; utility shape can favour balanced compromises over extremes. |
| 6 | Goal Programming - classical | Minimise weighted positive deviations from target values. | An unreachable goal remains missed; explain why raw and normalised deviation weights differ. |
| 7 | Lexicographic - classical | Optimise objectives in priority order and constrain later stages to preserve earlier optima. | Priority order is decisive; a lower-priority objective cannot trade away a higher-priority result. |
| 8 | TOPSIS - MCDM ranking | Rank alternatives by relative distance to a positive ideal and a negative ideal. | It selects from a given alternative set; criterion weights can change the order. |
| 9 | VIKOR - MCDM ranking | Rank compromise solutions using group utility and worst-criterion regret. | Parameter v controls the utility/regret balance; this is a ranker, not a Pareto optimizer. |
| 10 | VEGA - evolutionary | Form objective-specific mating subpopulations, then recombine their members. | Selection can over-represent extremes and under-represent balanced trade-offs. |
| 11 | MOGA - evolutionary | Rank by dominance and use fitness sharing to maintain diversity. | The project version is non-elitist; good solutions can be lost between generations. |
| 12 | NSGA - evolutionary | Use non-dominated ranks and sharing/crowding ideas to select a population. | Distinguish this earlier scheme from NSGA-II's explicit elitism and crowding selection. |
| 13 | NPGA - evolutionary | Compare candidates through Pareto-dominance tournaments and niche pressure. | Explain how tournament sampling and niche rules affect diversity and selection pressure. |
| 14 | NSGA-II - elitist evolutionary | Combine parents and offspring; select by non-dominated rank, then crowding distance. | This catalogue entry deliberately disables feasibility rules; compare with ID 21. |
| 15 | NSGA-III - many-objective evolutionary | Associate candidates with reference directions to preserve coverage as objective count rises. | Reference-point placement and objective scaling affect niche occupancy. |
| 16 | SPEA2 - elitist evolutionary | Combine dominance strength, raw fitness, density, and an external archive. | Archive truncation and density estimation influence diversity; runtime can exceed simpler methods. |
| 17 | PESA-II - elitist evolutionary | Select from a bounded archive using an objective-space hyper-grid, favouring sparse boxes. | Grid resolution controls the diversity-versus-selection balance. |
| 18 | MOEA/D - decomposition | Optimise a set of scalar subproblems whose weight vectors cover the objective simplex. | Neighbourhood information couples subproblems; weight coverage determines front coverage. |
| 19 | MOEA/D-DE - decomposition | Use MOEA/D subproblems with differential-evolution variation. | Performance depends on decomposition weights and DE variation parameters. |
| 20 | MOEA/D-ACO - discrete decomposition | Give each scalar subproblem pheromone-guided tour construction and update neighbouring incumbents. | Evaporation rho trades trail memory against adaptation; this is a compact educational variant. |
| 21 | Constrained NSGA-II | Apply Deb-style feasibility-first comparison, then Pareto rank and crowding. | Feasible candidates beat infeasible ones; among infeasible candidates, lower violation wins. |
| 22 | MOPSO - swarm | Move particles using velocity, personal best, and leaders from a non-dominated archive. | Archive leader choice is central to diversity; velocity and boundary handling affect stability. |
| 23 | MOGWO - swarm | Move candidates using three archive leaders analogous to alpha, beta, and delta wolves. | Leader quality and the shrinking exploration parameter control convergence and spread. |
| 24 | MOABC - swarm | Use employed, onlooker, and scout phases with dominance-aware greedy acceptance. | The trial limit triggers scouts; explain how this restores exploration after stagnation. |
| 25 | MOACO - discrete swarm | Construct tours with preference-weighted pheromone and heuristic information for objectives. | Separate pheromone trails encode competing criteria; archive quality drives reinforcement. |
| 26 | MOSA - physics-inspired | Use simulated-annealing chains with an archive and temperature-based acceptance. | Cooling governs exploration; the project uses a compact multi-chain archive version. |
| 27 | MODE - evolutionary | Apply DE/rand/1/bin and replace a target when its trial dominates it. | Dominance-only replacement can stall when many candidates are mutually non-dominated. |
| 28 | MO-TLBO - human-inspired | Alternate teacher and learner phases around a selected non-dominated teacher. | Teaching factor and acceptance rules shape exploration; this is a compact variant. |
| 29 | MOBA - swarm | Use frequency-driven motion, archive leaders, local walks, loudness, and pulse rate. | Archive leader selection affects coverage; the project uses a compact bat-inspired update. |
| 30 | MOMA - hybrid | Combine NSGA-II global search with short local hill-climbs under random scalar weights. | Local search improves selected candidates but consumes evaluations and can narrow preferences. |
| 31 | RVEA - many-objective | Associate candidates with reference vectors and use angle-penalised distance for selection. | Reference-vector adaptation and the penalty schedule matter as objective count grows. |
| 32 | IBEA - indicator | Assign fitness from pairwise additive epsilon-indicator values, then iteratively remove the worst. | Kappa controls indicator sensitivity; report both HV and IGD when tuning it. |
| 33 | SMS-EMOA - indicator | Use steady-state replacement, removing a candidate with low hypervolume contribution. | HV contribution preserves useful front coverage but can be computationally expensive. |
| 34 | HypE - indicator | Estimate hypervolume contribution with Monte Carlo samples for selection. | More samples reduce estimation noise at higher cost; multi-objective HV here is approximate. |
| 35 | R-NSGA-II - preference | Use NSGA-II ranking, then prefer candidates near a decision-maker reference point. | The result is preference-focused, not a uniform approximation of the full front. |
| 36 | MGDA - gradient | Find a minimum-norm convex combination of task gradients and descend along the common direction. | At a Pareto-stationary point the common descent direction vanishes; loss scaling changes behavior. |
| 37 | GradNorm - gradient | Adapt task weights so gradient magnitudes track relative task training rates. | Weight adaptation can oscillate and depends on step sizes and task-loss scales. |

## Quick comparisons

- **IDs 1 and 4:** both sweep scalar preferences; the weighted sum uses a linear combination, while Tchebycheff minimises the largest weighted deviation.
- **IDs 14 and 21:** same NSGA-II machinery in this code, but ID 14 ignores feasibility and ID 21 applies Deb-style constraint handling.
- **IDs 18-20:** MOEA/D decomposes the objective space; ID 19 changes the variation operator to DE, while ID 20 couples decomposition to ant-colony tour construction.
- **IDs 32-34:** indicator-based selection uses epsilon, hypervolume contribution, or an estimated hypervolume contribution, respectively.
- **IDs 8 and 9:** decision-support ranking methods; they rank a finite alternative set and do not generate an approximation set.
- **IDs 36 and 37:** gradient methods for differentiable P13(a), not general black-box optimisers in this implementation.

## Evidence links

- Experiment summaries: `results/implementation_log.csv` and `results/taskD_compare.txt`.
- Task C trace and sensitivity: `results/handtrace_taskC.txt` and `results/taskC_sensitivity.txt`.
- Report with tables, interpretations, and figure captions: `REPORT.html` and `REPORT.pdf`.