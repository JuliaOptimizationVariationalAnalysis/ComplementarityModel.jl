using ComplementarityModel
using Documenter

DocMeta.setdocmeta!(
  ComplementarityModel,
  :DocTestSetup,
  :(using ComplementarityModel);
  recursive = true,
)

makedocs(;
  modules = [ComplementarityModel],
  authors = "Tangi Migot tangi.migot@gmail.com",
  repo = "https://github.com/JuliaOptimizationVariationalAnalysis/ComplementarityModel.jl/blob/{commit}{path}#{line}",
  sitename = "ComplementarityModel.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://JuliaOptimizationVariationalAnalysis.github.io/ComplementarityModel.jl",
    assets = String[],
  ),
  pages = ["Home" => "index.md"],
)

deploydocs(;
  repo = "github.com/JuliaOptimizationVariationalAnalysis/ComplementarityModel.jl",
  devbranch = "main",
)
