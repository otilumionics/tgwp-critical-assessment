module Spectra
using FFTW
using ..GaussianToolbox

function autocorrelation(t,
                         p0::Real, q0::Real, u0::Real, w0::Real, λ0::Real,
                         p::Real,  q::Real,  u::Real,  w::Real,  λ::Real)::Complex
    N0 = sqrt(gaussian_norm2(p0, q0, u0, w0, λ0))
    Nt = sqrt(gaussian_norm2(p, q, u, w, λ))

    return gaussian_overlap(p0, q0, u0, w0, λ0, p,  q,  u,  w,  λ)/(N0*Nt)
end

@inline spectrum(C) = fftshift(ifft(C))

@inline freqs(t_arr, dt; scale=1, shift=0) =
    scale*fftshift(fftfreq(length(t_arr), 1/dt)) .+ shift

end # module Spectra
