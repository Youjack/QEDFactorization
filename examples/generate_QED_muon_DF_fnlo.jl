using QEDFactorization.QEDEvolve
using Printf

const setname = "QED_muon_DF_fnlo"
const id = 13

const me = 0.000511
const mμ = 0.105658
const αEMe = 1 / 137.036
const αEMμ = evolve_αEM_from_μ₀²_to_μ²(αEMe, me^2, mμ^2, 1)

# τ is not considered
function get_αEM(μ²::Real)::Float64
    evolve_αEM_from_μ₀²_to_μ²(αEMμ, mμ^2, μ², 2)
end

# Q grid copied from JAM22 (same as JAM20)
const Q_grid = [ 1.14018E+00, 1.22539E+00, 1.32482E+00, 1.44156E+00, 1.57954E+00, 1.74381E+00, 
                 1.94091E+00, 2.12132E+00, 2.17945E+00, 2.49622E+00, 2.89383E+00, 3.39925E+00, 
                 4.05063E+00, 4.90286E+00, 6.03623E+00, 7.57063E+00, 9.68869E+00, 1.26749E+01, 
                 1.69835E+01, 2.33578E+01, 3.30501E+01, 4.82334E+01, 7.28041E+01, 1.13998E+02, 
                 1.85779E+02, 3.16228E+02
]
# x grid copied from JAM22 (same as JAM20)
const x_grid = [ 1.00000E-06, 1.28121E-06, 1.64152E-06, 2.10317E-06, 2.69463E-06, 3.45242E-06, 
                 4.42329E-06, 5.66715E-06, 7.26076E-06, 9.30241E-06, 1.19180E-05, 1.52689E-05, 
                 1.95617E-05, 2.50609E-05, 3.21053E-05, 4.11287E-05, 5.26863E-05, 6.74889E-05, 
                 8.64459E-05, 1.10720E-04, 1.41800E-04, 1.81585E-04, 2.32503E-04, 2.97652E-04, 
                 3.80981E-04, 4.87518E-04, 6.26039E-04, 8.00452E-04, 1.02297E-03, 1.30657E-03, 
                 1.66759E-03, 2.12729E-03, 2.71054E-03, 3.44865E-03, 4.37927E-03, 5.54908E-03, 
                 7.01192E-03, 8.83064E-03, 1.10763E-02, 1.38266E-02, 1.71641E-02, 2.11717E-02, 
                 2.59364E-02, 3.15062E-02, 3.79623E-02, 4.53425E-02, 5.36750E-02, 6.29705E-02, 
                 7.32221E-02, 8.44039E-02, 9.64793E-02, 1.09332E-01, 1.23067E-01, 1.37507E-01, 
                 1.52639E-01, 1.68416E-01, 1.84794E-01, 2.01731E-01, 2.19016E-01, 2.36948E-01, 
                 2.55242E-01, 2.73927E-01, 2.92954E-01, 3.12340E-01, 3.32036E-01, 3.52019E-01, 
                 3.72282E-01, 3.92772E-01, 4.13533E-01, 4.34326E-01, 4.55495E-01, 4.76836E-01, 
                 4.98342E-01, 5.20006E-01, 5.41818E-01, 5.63773E-01, 5.85861E-01, 6.08077E-01, 
                 6.30459E-01, 6.52800E-01, 6.75387E-01, 6.98063E-01, 7.20830E-01, 7.43683E-01, 
                 7.66623E-01, 7.89636E-01, 8.12791E-01, 8.35940E-01, 8.59175E-01, 8.82485E-01, 
                 9.05866E-01, 9.29311E-01, 9.52817E-01, 9.76387E-01, 9.99999E-01
]

function main()

info_str = """
SetDesc:         "QED distribution functions of muon, fixed at NLO."
Authors:         Lin Youjack
Reference:       N/A
Format:          lhagrid1
DataVersion:     1
NumMembers:      1
SetIndex:        00000
Particle:        $id
Flavors:         [$id, $(100id)]
OrderQCD:        1
NumFlavors:      3
XMin:            1.00E-06
XMax:            1.00E+00
QMin:            1.14E+00
QMax:            3.16E+02
AlphaS_Type:     ipol
AlphaS_OrderQCD: 1
"""

# Q, x grid ----------------------------------------------------------------------------------------

info_str *= "AlphaS_Qs: [ "
for Q in Q_grid
    info_str *= @sprintf("%.5E, ", Q)
end
info_str = info_str[1:end-2]
info_str *= " ]\n"
info_str *= "AlphaS_Vals: [ "
for Q in Q_grid
    info_str *= @sprintf("%.5E, ", get_αEM(Q^2))
end
info_str = info_str[1:end-2] * " ]\n"

dat_head_str = """
PdfType: central
Format: lhagrid1
---
"""
for x in x_grid
    dat_head_str *= @sprintf("%.5E ", x)
end
dat_head_str = dat_head_str[1:end-1] * "\n"
for Q in Q_grid
    dat_head_str *= @sprintf("%.5E ", Q)
end
dat_head_str = dat_head_str[1:end-1] * "\n"
dat_head_str *= "$id $(100id)\n"

# distribution functions ---------------------------------------------------------------------------

# μ > mμ
f̂eμ²(μ²) = s -> let αEM = get_αEM(μ²)
    f̂ll₀(αEM,s) + αEM/2π * QEDEvolve.P̂ll(s) * log(μ²/mμ^2)
end

# write files --------------------------------------------------------------------------------------

folder = joinpath(split(ENV["LHAPDF_DATA_PATH"], ":")[1], setname)
if !isdir(folder) mkdir(folder) end

open(joinpath(folder, "$(setname).info"), "w") do info_file
    write(info_file, info_str)
end

open(joinpath(folder, "$(setname)_0000.dat"), "w") do dat_file
    write(dat_file, dat_head_str)

    for x in x_grid
        for Q in Q_grid
            write(dat_file,
                @sprintf("%.5E", (1-x) * mellin_inv(f̂eμ²(Q^2))(1-x)         ) * " " *
                @sprintf("%.5E", 1 - integrate_distribution(f̂eμ²(Q^2), 1-x) ) * "\n" )
        end
    end
    
    write(dat_file, "---\n")
end

end # main()

main()
