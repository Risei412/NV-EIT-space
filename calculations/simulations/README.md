# calculations/simulations/

Reserved for dynamics and Monte Carlo campaigns.

Currently empty. The two campaigns that might have landed here did not:

- the **room-temperature campaign** (`RoomT/`, steps 1–9) is a parameter-sweep
  and global-optimization campaign over a static Liouvillian, so it stayed in
  `../numerics/New no-go theory/RoomT/`;
- the **Monte Carlo confidence intervals** reported by Gate E are computed
  inside that gate's runner rather than as a standalone simulation.

Put a campaign here when it integrates dynamics in time or samples
trajectories.
