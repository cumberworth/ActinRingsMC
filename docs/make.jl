using ActinRingsMC
using Documenter

DocMeta.setdocmeta!(ActinRingsMC, :DocTestSetup, :(using ActinRingsMC); recursive=true)

makedocs(;
    modules=[ActinRingsMC],
    authors="Alexander Cumberworth <alex@cumberworth.org>",
    repo="https://github.com/cumberworth/ActinRingsMC.jl/blob/{commit}{path}#{line}",
    sitename="ActinRingsMC.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cumberworth.github.io/ActinRingsMC.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cumberworth/ActinRingsMC.jl",
    devbranch="master",
)
