using argmin
using Documenter

makedocs(;
    modules=[argmin],
    authors="Atsushi Sakai <asakai.amsl+github@gmail.com> and contributors",
    repo="https://github.com/AtsushiSakai/argmin.jl/blob/{commit}{path}#L{line}",
    sitename="argmin.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://AtsushiSakai.github.io/argmin.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AtsushiSakai/argmin.jl",
)
