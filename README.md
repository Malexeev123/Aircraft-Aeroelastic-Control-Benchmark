# Aircraft-Aeroelastic-Control-Benchmark

**This project is  _under construction_** 

This comprehensive benchmark allows to couple a MATLAB
control framework to the SHARPy nonlinear aeroelastic
solver, effectively creating a virtual hardware-in-the-loop environment.

## MATLAB Repository structure:
###   /rom: structural and aerodynamic ROM data; modal maps; gravity and follower-force mappings.
###    /plant: time-stepper, outputs, sensor models.
###    /est: MHE problem; arrival covariance; window management; warm starts.
###    /ctrl: MPC multiple shooting; constraints; warm starts; LQR design.
###    /fusion: complementary filters and interface to actuator chain.
###    /act: servo model, saturation, rate limits, anti-windup.
 ###   /cases: A/B/C scripts with identical seeds and logging.
 ###   /post: metrics (GLA\%, RMS, DEL), actuator statistics, timing tables, plots.
 ###   /validate: high-fidelity cross-check scripts and plotting templates.

