#files = [
#		]

#for f in files
#    include(f * ".jl")
#end

module constants
    using Unitful
    gEarth = 9.80665u"m/s^2"
    GEarth = (6.67408*10^(-11))u"N*m^2/kg^2"
    k = (8.9875517887*10^9)u"N*m^2/C^2"
    c = 299792458u"m/s"
    e = (1.60217662*10^(-19))u"C"
    h = (6.62607004*10^(-34))u"J*s"
    eps0 = (8.85418782*10^(-12))u"m^(-3)*kg^(-1)*s^(4)*A^(2)"
    mu0 = (1.25663706*10^(-6))u"m*kg*s^-2*A^-2"
    mEarth = (5.972*10^24)u"kg"
    rEarth = 6378.1u"km"
    mSun = (1.989*10^30)u"kg"
    export
        gEarth,
        GEarth,
        k,
        c,
        e,
        h,
        eps0,
        mu0,
        mEarth,
        rEarth,
        mSun
end

module physics
    module kinematics
        using ....constants
        # Displacement
        x_abs(x0,x) = abs(x-x0)
        x_const_acc1(v0,a,t) = v0*t+1/2*a*t^2
        x_const_acc2(v0,v,t) = ((v+v0)/2)*t
        # Velocity
        v(x0,x,t0,t) = (x-x0)/(t-t0)
        v_const_acc1(v0,a,t) = v0+a*t
        v_const_acc2(x0,x,v0,a) = sqrt(v0^2+2*a*(x-x0))
        # Acceleration
        a(x0,x,v0,v) = (v^2-v0^2)/(2*(x-x0))
        R(v0,θ0) = v0^2/gEarth*sin(2*θ0)
        y_max(v0,θ0) = v0^2/(2*gEarth)*sin(θ0)^2
    end

    module dynamics
        using ....constants
        using LinearAlgebra
        # Newton's Laws
        F(m,a) = m*a
        F_reaction(F_in) = -F_in
        # Forces
        F_weight(m,g) = m*g
        F_friction(u,N) = u*N
        F_centripetal(m,v,r) = m*v^2/r
        # Work and Energy
        W(F,d) = dot(F,d)
        W_KE(m,v0,v) = 1/2*m*(v^2-v0^2)
        U_g(m,g,h) = m*g*h
        P_avg(W0,W,t0,t) = (W-W0)/(t-t0)
        P_inst(F,v) = dot(F,v)
        p_linear(m,v) = m*v
        I(F,t0,t) = F*(t-t0)
        # Gravitation
        F_g(G,m1,m2,r) = G*m1*m2/r^2
        g(G,M,r) = G*M/r^2
        U(G,m1,m2,r) = -G*m1*m2/r
        v_escape(G,m,r) = sqrt(2*G*m/r)
        T(a,G,M) = sqrt(4*pi^2*a^3/(G*M))
        a(T,G,M) = (G*M*T^2/(4*pi^2))^(1/3)
    end

    module rotationalDynamics
        using LinearAlgebra
        # Uniform Circular Motion
        θ(s,r) = s/r
        ω(v_t,r) = v_t/r
        # Period
        T_f(f) = 1/f
        T_ω(ω) = 2*pi/ω
        # Average Angular Acceleration
        α_avg(a_t,r) = a_t/r
        # Centripetal Accleration
        a_c(v,t,r) = v*t^2/r
        # Rotational Motion with Constant Acceleration
        ω_const_acc1(ω_0,α,t) = ω_0+α*t
        ω_const_acc2(ω_0,α,θ_0,θ) = sqrt(ω_0^2+2*α*(θ-θ_0))
        θ_const_acc1(θ_0,ω_0,α,t) = θ_0+ω_0*t+1/2*α*t^2
        θ_const_acc2(θ_0,ω_0,ω,t) = θ_0+1/2*(ω_0+ω)*t
        # Torque
        τ(r,F) = cross(r,F)
        τ_net(I,α) = I*α
        # Angular Momentum
        L(r,p) = cross(r,p)
        L_mag(r,m,v,θ) = r*m*v*sin(θ)
        L_rigid(I,ω) = I*ω
        # Moment of Inertia
        function I(m_arr,r_arr)
            I = 0
            if length(m_arr) != length(r_arr)
                return "ERR"
            end
            for i = 1:length(m_arr)
                I+=m_arr[i]*r_arr[i]^2
            end
            return I
        end
        function I_common(m,R,l,body_type)
            x = body_type
            if x == "particle" || x == "cylinder-thin" || x == "hoop-thin"
                I = m*R^2
            elseif x == "cylinder-xy" || x == "disc-xy" || x == "hoop-xy"
                I = 1/2*m*R^2
            elseif x == "disc-thin"
                I = 1/4*m*R^2
            elseif x == "sphere-solid"
                I = 2/5*m*R^2
            elseif x == "sphere-hollow"
                I = 2/3*m*R^2
            elseif x == "rod-center"
                I = 1/12*m*l^2
            elseif x == "rod-end"
                I = 1/3*m*l^2
            end
        end
        K(I,ω) = 1/2*I*ω^2
        W(τ,θ) = τ*θ
        P(τ,ω) = τ*ω
    end

    module simpleHarmonics
        # Spring with a Mass
        F(k,x) = -k*x
        f(k,m) = 1/(2*pi)*sqrt(k/m)
        ω(k,m) = sqrt(k/m)
        T(k,m) = (2*pi)*sqrt(m/k)
        # Equations of Motion
        # A is amplitude of oscillation
        x(A,ω,t) = A*cos(ω*t)
        v(A,ω,t) = -ω*A*sin(ω*t)
        a(A,ω,t) = -ω^2*A*cos(ω*t)
        # Energy Equations
        U(k,x) = 1/2*k*x^2
        E(k,A) = 1/2*k*A^2
        # Pendulums
        T_pend(l,g) = 2*pi*sqrt(l/g)
    end
end