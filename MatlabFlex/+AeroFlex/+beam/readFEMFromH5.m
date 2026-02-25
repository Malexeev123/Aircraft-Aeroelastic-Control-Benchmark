function fem = readFEMFromH5(h5file)
    fem.num_node_elem = h5read(h5file,'/num_node_elem');
    fem.num_elem      = h5read(h5file,'/num_elem');
    fem.num_node      = h5read(h5file,'/num_node');
    fem.coordinates   = double(h5read(h5file,'/coordinates')); 
    fem.connectivities= double(h5read(h5file,'/connectivities'));
    fem.stiffness_db  = double(h5read(h5file,'/stiffness_db'));
    fem.elem_stiffness= double(h5read(h5file,'/elem_stiffness'));
    fem.mass_db       = double(h5read(h5file,'/mass_db'));
    fem.elem_mass     = double(h5read(h5file,'/elem_mass'));
    fem.boundary_conditions = double(h5read(h5file,'/boundary_conditions'));
    fem.app_forces = double(h5read(h5file,'/app_forces'));
    fem.beam_number = double(h5read(h5file,'/beam_number'));
    if h5exists(h5file, '/lumped_mass_nodes')
        fem.lumped_mass = double(h5read(h5file,'/lumped_mass'));
        fem.lumped_mass_nodes = h5read(h5file,'/lumped_mass_nodes');
        fem.lumped_mass_inertia = double(h5read(h5file,'/lumped_mass_inertia'));
        fem.lumped_mass_position= double(h5read(h5file,'/lumped_mass_position'));
    else
        fem.lumped_mass = [];
        fem.lumped_mass_nodes = [];
        fem.lumped_mass_inertia=[];
        fem.lumped_mass_position=[];
    end
    szPos = size(fem.lumped_mass_position);
    if szPos(1) == 3 && szPos(2) == length(fem.lumped_mass)
        fem.lumped_mass_position = fem.lumped_mass_position.'; % => [nL,3]
    end
    
    if h5exists(h5file, '/structural_twist')
       fem.structural_twist = double(h5read(h5file,'/structural_twist'));
    else
       fem.structural_twist = zeros(fem.num_elem,1);
    end
    if h5exists(h5file,'/frame_of_reference_delta')
       fem.frame_of_reference_delta= double(h5read(h5file,'/frame_of_reference_delta'));
    else
       fem.frame_of_reference_delta= zeros(fem.num_elem,3);
    end
    
    [rC,cC] = size(fem.coordinates);
    if rC==3 && cC==fem.num_node
       fem.coordinates= fem.coordinates';
    end
    [rE,cE]= size(fem.connectivities);
    if rE~=fem.num_elem && cE==fem.num_elem
       fem.connectivities= fem.connectivities';
    end
    fem.connectivities= fem.connectivities+1;
    end
    
    function flag = h5exists(h5file, ds)
    info = h5info(h5file);
    allp = recList(info,'');
    flag= ismember(ds, allp);
    end
    function lst = recList(info, prefix)
    lst= {};
    for i=1:length(info.Datasets)
       nm= [prefix '/' info.Datasets(i).Name];
       lst{end+1}= nm; 
    end
    for ig=1:length(info.Groups)
       sp= [prefix '/' info.Groups(ig).Name];
       subn= recList(info.Groups(ig), sp);
       lst= [lst, subn]; 
    end
end