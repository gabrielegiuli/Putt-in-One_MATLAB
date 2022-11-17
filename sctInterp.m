function z_interp = sctInterp(x_sct, y_sct, z_sct, x0, y0)
    z_interp = griddata(x_sct, y_sct, z_sct, x0, y0, 'linear');
end