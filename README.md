# ReputationAwareUninormDriven
Reputation-based blockchain consensus framework using Intuitionistic Fuzzy Sets (IFS) and Uninorm Aggregation Operators (UAO). The model computes a dynamic reputation degree and reputation weight to manage validator behaviour under uncertainty, enabling recovery from occasional errors while penalising repeated misbehaviour. The strictness of evaluation is controlled via a tunable parameter (α). The approach maintains linear computational complexity and introduces no additional communication overhead. The framework models validator behaviour under uncertainty through a two-phase process:
## Reputation Degree (RepD) Computation
- Derived from the Successful Validation Rate (SuccVR).
- Uses intuitionistic fuzzy components: membership (μ), non-membership (ν), and hesitation (π).
- Produces a bounded and monotonic reputation score reflecting short-term validator performance.
## Reputation Weight (w) Aggregation
- Applies a uninorm operator (e.g., Fodor’s uninorm) to aggregate current reputation degree with historical reputation weight.
- Captures long-term trust dynamics.
- Reinforces consistent behaviour and penalises repeated misbehaviour.
## Key Features
- Models uncertainty explicitly using IFS theory.
- Supports reputation recovery after isolated errors.
- Penalises consecutive faults through uninorm-based reinforcement.
- Tunable strictness via parameter α, controlling sensitivity of reputation evaluation.
- Linear computational complexity per round: O(n) for n validators.
- No additional communication overhead beyond the underlying consensus protocol.
