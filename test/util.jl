
function ensure_vector(n, ue::Vector, ui::Vector, ua::Vector)
    return ue, ui, ua
end

function ensure_vector(n, ue::Real, ui::Real, ua::Real)
    ue = ones(n)*ue
    ui = ones(n)*ui
    ua = ones(n)*ua
    return ue, ui, ua
end

function dx_step(p, prev_ue, prev_ui, prev_ua)
    n = p.n

    prev_ue, prev_ui, prev_ua = ensure_vector(n, prev_ue, prev_ui, prev_ua)

    initials = zeros(3*n)
    initials[1:n] .= prev_ue
    initials[n+1:2*n] .= prev_ui
    initials[2*n+1:end] .= prev_ua
    hist(x, p, t) = (x .= initials)
    
    dx = zeros(3*n)
    WCL.dx!(dx, initials, hist, p, 2.0)
    return dx
end

function prepare_euler_arr(n, ue, ui, ua)
    ue, ui, ua = ensure_vector(n, ue, ui, ua)

    ue_arr = zeros(n, 5)
    ui_arr = zeros(n, 5)
    ua_arr = zeros(n, 5)

    for i in 1:5
        ue_arr[:, i] .= ue
        ui_arr[:, i] .= ui
        ua_arr[:, i] .= ua
    end

    return ue_arr, ui_arr, ua_arr
end

function dx_euler_step(p, prev_ue, prev_ui, prev_ua)
    ue, ui, ua = prepare_euler_arr(p.n, prev_ue, prev_ui, prev_ua)
    
    due = zeros(p.n)
    dui = zeros(p.n)
    dua = zeros(p.n)

    WCL.dx_euler!(due, dui, dua, ue, ui, ua, 0.1, p, 2)

    return due, dui, dua
end

function dx_form_euler_step(p, prev_ue, prev_ui, prev_ua)
    ue, ui, ua = prepare_euler_arr(p.n, prev_ue, prev_ui, prev_ua)
    
    due = zeros(p.n)
    dui = zeros(p.n)
    dua = zeros(p.n)

    WCL.dx_matrix_euler!(due, dui, dua, ue, ui, ua, p, 2)

    return due, dui, dua
end