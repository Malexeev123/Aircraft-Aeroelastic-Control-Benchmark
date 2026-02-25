function aero = readAEROFromH5(h5file)
    aero.aero_node = h5read(h5file,'/aero_node');
    aero.airfoil_distribution      = h5read(h5file,'/airfoil_distribution');
    aero.chord   = (h5read(h5file,'/chord')); 
    aero.control_surface= (h5read(h5file,'/control_surface'));
    aero.control_surface_chord  = (h5read(h5file,'/control_surface_chord'));
    aero.control_surface_deflection= (h5read(h5file,'/control_surface_deflection'));
    aero.control_surface_hinge_coord       = (h5read(h5file,'/control_surface_hinge_coord'));
    aero.control_surface_type     = (h5read(h5file,'/control_surface_type'));
    aero.elastic_axis = (h5read(h5file,'/elastic_axis'));
    aero.m_distribution = (h5read(h5file,'/m_distribution'));
    aero.surface_distribution = (h5read(h5file,'/surface_distribution'));
    aero.surface_m = (h5read(h5file,'/surface_m'));
    aero.twist = (h5read(h5file,'/twist'));
    
    
    
    aero.airfoils.one      = h5read(h5file,'/airfoils/0');
    aero.airfoils.two      = h5read(h5file,'/airfoils/1');
    aero.airfoils.three      = h5read(h5file,'/airfoils/2');
    aero.airfoils.four      = h5read(h5file,'/airfoils/3');
    aero.airfoils.five      = h5read(h5file,'/airfoils/4');
    aero.airfoils.six      = h5read(h5file,'/airfoils/5');
    aero.airfoils.seven      = h5read(h5file,'/airfoils/6');
    aero.airfoils.eight      = h5read(h5file,'/airfoils/7');
    aero.airfoils.nine      = h5read(h5file,'/airfoils/8');
end