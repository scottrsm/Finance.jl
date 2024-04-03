using Finance
import Pkg

Pkg.add("Documenter")
using Documenter

makedocs(
	sitename = "Finance",
	format = Documenter.HTML(),
	modules = [Finance]
	)

	# Documenter can also automatically deploy documentation to gh-pages.
	# See "Hosting Documentation" and deploydocs() in the Documenter manual
	# for more information.
	deploydocs(
		repo = "github.com/scottrsm/Finance.jl.git"
	)
