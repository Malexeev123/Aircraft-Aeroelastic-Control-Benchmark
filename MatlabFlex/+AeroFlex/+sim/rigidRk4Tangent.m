function [next,stateTangent,wrenchTangent] = ...
        rigidRk4Tangent(state,wrench,parameters,dt)
%RIGIDRK4TANGENT Advance Newton-Euler dynamics and exact RK4 tangents.
%   The state order is [r_I; v_B; euler_321; omega_B]. The input order is
%   [F_B; M_B], with force in N and moment in N m.

state = state(:);
wrench = wrench(:);
assert(numel(state)==12 && numel(wrench)==6 && ...
    isstruct(parameters) && isfield(parameters,"mass") && ...
    isfield(parameters,"I_B") && isscalar(parameters.mass) && ...
    isfinite(parameters.mass) && parameters.mass>0 && ...
    isequal(size(parameters.I_B),[3,3]) && ...
    all(isfinite(parameters.I_B),"all") && ...
    isscalar(dt) && isfinite(dt) && dt>0, ...
    "AeroFlex:RigidRk4TangentInput", ...
    "Rigid RK4 tangent inputs are incomplete or invalid.");
assert(all(isfinite(real(state))) && all(isfinite(imag(state))) && ...
    all(isfinite(real(wrench))) && all(isfinite(imag(wrench))), ...
    "AeroFlex:RigidRk4TangentFinite", ...
    "Rigid RK4 tangent inputs must be finite.");

identity = eye(12,"like",state+wrench(1));
zeroInput = zeros(12,6,"like",identity);

[k1,A1,B1] = localDerivative(state,wrench,parameters);
stateStage2 = state+0.5*dt*k1;
stateSensitivity2 = identity+0.5*dt*A1;
inputSensitivity2 = zeroInput+0.5*dt*B1;

[k2,A2,B2] = localDerivative(stateStage2,wrench,parameters);
k2State = A2*stateSensitivity2;
k2Input = A2*inputSensitivity2+B2;
stateStage3 = state+0.5*dt*k2;
stateSensitivity3 = identity+0.5*dt*k2State;
inputSensitivity3 = zeroInput+0.5*dt*k2Input;

[k3,A3,B3] = localDerivative(stateStage3,wrench,parameters);
k3State = A3*stateSensitivity3;
k3Input = A3*inputSensitivity3+B3;
stateStage4 = state+dt*k3;
stateSensitivity4 = identity+dt*k3State;
inputSensitivity4 = zeroInput+dt*k3Input;

[k4,A4,B4] = localDerivative(stateStage4,wrench,parameters);
k4State = A4*stateSensitivity4;
k4Input = A4*inputSensitivity4+B4;

next = state+(dt/6)*(k1+2*k2+2*k3+k4);
stateTangent = identity+(dt/6)*(A1+2*k2State+2*k3State+k4State);
wrenchTangent = (dt/6)*(B1+2*k2Input+2*k3Input+k4Input);
end

function [value,stateJacobian,inputJacobian] = ...
        localDerivative(state,wrench,parameters)
vBody = state(4:6);
euler = state(7:9);
omegaBody = state(10:12);
forceBody = wrench(1:3);
momentBody = wrench(4:6);
inertia = parameters.I_B;

[rotation,rotationDerivative] = localBodyToInertial(euler);
[rateMap,rateMapDerivative] = localEulerRateMap(euler);
value = [rotation*vBody; ...
    forceBody/parameters.mass-localSkew(omegaBody)*vBody; ...
    rateMap*omegaBody; ...
    inertia\(momentBody-localSkew(omegaBody)*(inertia*omegaBody))];

stateJacobian = zeros(12,12,"like",value);
stateJacobian(1:3,4:6) = rotation;
for coordinate = 1:3
    stateJacobian(1:3,6+coordinate) = ...
        rotationDerivative(:,:,coordinate)*vBody;
end
stateJacobian(4:6,4:6) = -localSkew(omegaBody);
stateJacobian(4:6,10:12) = localSkew(vBody);
stateJacobian(7:9,10:12) = rateMap;
for coordinate = 1:2
    stateJacobian(7:9,6+coordinate) = ...
        rateMapDerivative(:,:,coordinate)*omegaBody;
end
stateJacobian(10:12,10:12) = inertia\( ...
    localSkew(inertia*omegaBody)-localSkew(omegaBody)*inertia);

inputJacobian = zeros(12,6,"like",value);
inputJacobian(4:6,1:3) = eye(3)/parameters.mass;
inputJacobian(10:12,4:6) = inertia\eye(3);
end

function [rotation,derivative] = localBodyToInertial(euler)
phi=euler(1); theta=euler(2); psi=euler(3);
cp=cos(phi); sp=sin(phi); ct=cos(theta); st=sin(theta);
cs=cos(psi); ss=sin(psi);
rotation=[ct*cs,sp*st*cs-cp*ss,cp*st*cs+sp*ss; ...
    ct*ss,sp*st*ss+cp*cs,cp*st*ss-sp*cs; ...
    -st,sp*ct,cp*ct];
derivative=zeros(3,3,3,"like",euler);
derivative(:,:,1)=[0,cp*st*cs+sp*ss,-sp*st*cs+cp*ss; ...
    0,cp*st*ss-sp*cs,-sp*st*ss-cp*cs;0,cp*ct,-sp*ct];
derivative(:,:,2)=[-st*cs,sp*ct*cs,cp*ct*cs; ...
    -st*ss,sp*ct*ss,cp*ct*ss;-ct,-sp*st,-cp*st];
derivative(:,:,3)=[-ct*ss,-sp*st*ss-cp*cs,-cp*st*ss+sp*cs; ...
    ct*cs,sp*st*cs-cp*ss,cp*st*cs+sp*ss;0,0,0];
end

function [map,derivative] = localEulerRateMap(euler)
phi=euler(1); theta=euler(2);
sp=sin(phi); cp=cos(phi); tangent=tan(theta); ct=cos(theta);
assert(abs(ct)>1e-8,"AeroFlex:RigidRk4TangentEulerSingularity", ...
    "The analytical rigid tangent is undefined at the Euler singularity.");
map=[1,sp*tangent,cp*tangent;0,cp,-sp;0,sp/ct,cp/ct];
derivative=zeros(3,3,2,"like",euler);
derivative(:,:,1)=[0,cp*tangent,-sp*tangent;0,-sp,-cp; ...
    0,cp/ct,-sp/ct];
sec2=1/ct^2;
derivative(:,:,2)=[0,sp*sec2,cp*sec2;0,0,0; ...
    0,sp*sec2*sin(theta),cp*sec2*sin(theta)];
end

function matrix = localSkew(vector)
matrix = [0,-vector(3),vector(2);vector(3),0,-vector(1); ...
    -vector(2),vector(1),0];
end
