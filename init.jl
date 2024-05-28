isnothing(x) = x === nothing

section(x::String; bold=true, color=:default) =
    printstyled("%% ", x,'\n', bold=bold, color=color)

function mkfigure(k::Int=0)
    return function figure(;reset=false)
        return k = reset ? 1 : k+1
    end
end

macro echo(ex)
    return quote
        printstyled($(string(ex))*'\n', color=:light_black)
        local ans = $(esc(ex))
        print(" => ")
        show(stdout, "text/plain", ans)
        println()
        println()
        ans
    end
end
