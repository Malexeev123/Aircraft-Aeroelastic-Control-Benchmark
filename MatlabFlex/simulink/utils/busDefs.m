function busDefs
% Run once before opening any model
elems(1)=Simulink.BusElement('F'); elems(1).Dimensions=3;
elems(2)=Simulink.BusElement('M'); elems(2).Dimensions=3;
Simulink.Bus.addElementToBus(Simulink.Bus('ForceMomBus'),elems,true);

elems2(1)=Simulink.BusElement('q'); elems2(1).Dimensions=4;
elems2(2)=Simulink.BusElement('uvw'); elems2(2).Dimensions=3;
elems2(3)=Simulink.BusElement('pqr'); elems2(3).Dimensions=3;
Simulink.Bus.addElementToBus(Simulink.Bus('RBStateBus'),elems2,true);
end
