
## Meeting European conservation and restoration targets under future land-use demands

The European Union is committed to achieving ambitious area-based conservation and restoration targets in the upcoming decade. Yet, there is a contention that these targets may conflict with societal needs, particularly food and timber production. Ensuring that competing demands for land are met, while the underlying objectives to mitigate climate change and reduce biodiversity loss are achieved, will require strategic planning across land uses, management measures, and jurisdictional boundaries. Here, we use an integrated spatial planning approach to identify restoration and conservation measures that maximize benefits to species as well as contributions to climate mitigation, while ensuring projected 2030 crop, pasture, and forestry production demands are met across the EU.

### Code structure

All code used to format input data and run prioritization analysis is available in the `/scripts/problem-setup/` folder:

+ __`1-PU-globiom.R`:__ Sets up planning units and calculates initial land use proportions
+ __`2-zones.R` :__ Sets up 25 management zones and planning unit level constraints for each zone
+ __`3-features-updated.R` & `3b-feature-targets-disaggregated.R`:__ formats zone contributions to feature targets and sets disaggregated feature targets (by country-biome-spp)
+ __`4-linear-constraints-*.R`:__ sets constraints across nuts2 regions for f455 and BAU scenarios to ensure production targets are met at subnational scales ("linear constraints"). `
+ __`5-problem-solve-splitspp-*.R`:__ run spatial optimization problem for all scenarios. Scripts are split by future productions scenarios (Fit455 and Reference/BAU) and a separate script is used for the restoration-only scenarios (shown in Figure S3)
+ __`5b-format-solutions.R`:__ takes the outputs from the `problem-solve` scripts and calculates feature shortfalls, formats solutions as rasters, etc.
+ __`6-figures.R`:__ makes all figures panels used in manuscript

In this code, we rely on ![Gurobi v 9.5](https://www.gurobi.com/) to solve optimization solutions (a proprietary LP solver, offers free academic licenses). However, it is possible to solve these problems with alternative open-source solvers (more information on supported solvers and their benchmarks ![here](https://prioritizr.net/articles/solver_benchmarks.html)

### Data availability

All formatted data to create the optimization solutions presented in the paper (`/scripts/problem-setup/5-*.R`) is available in the `/data-formatted/` folder. Underlying data for setting up input files (steps 1-4) is available upon request.




