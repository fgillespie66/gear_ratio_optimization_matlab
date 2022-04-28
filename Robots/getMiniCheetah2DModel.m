function model = getMiniCheetah2DModel()
% This function builds a 2D version of the floating base model of the
% Mini Cheetah quadruped using the conventions of the spatialv2 library

%% Robot Parameters
bodyMass = 3.3;
bodyLength = 0.19 * 2;
bodyWidth = 0.049 * 2;
bodyHeight = 0.05 * 2;

abadLinkLength = 0.062;
hipLinkLength = 0.209;
kneeLinkLength = 0.195;
% kneeLinkY_offset = 0.004;

abadRotInertia = 1e-6 .* [381, 58, 0.45;...
    58, 560, 0.95;...
    0.45, 0.95, 444];
abadCOM = [0, 0.036, 0];
%abadInertia = spatialInertia(0.54, abadCOM, abadRotInertia);
abadInertia = mcI(0.54, abadCOM([1 3]), abadRotInertia(2,2) );

hipRotInertia = 1e-6 .* [1983, 245, 13;...
    245, 2103, 1.5;...
    13, 1.5, 408];
hipCOM = [0, 0.016, -0.02];
%hipInertia = spatialInertia(0.634, hipCOM, hipRotInertia);
hipInertia = mcI(0.634, hipCOM([1 3]), hipRotInertia(2,2) );

kneeRotInertia = 1e-6 .* [6, 0, 0;...
    0, 248, 0;...
    0, 0, 245];
kneeCOM = [0, 0, -0.061];
%kneeInertia = spatialInertia(0.064, kneeCOM, kneeRotInertia);
kneeInertia = mcI(0.064, kneeCOM([1 3]), kneeRotInertia(2,2) );

bodyRotInertia = 1e-6 .* [11253, 0, 0;...
    0, 36203, 0;...
    0, 0, 42673];
bodyCOM = [0, 0, 0];
%bodyInertia = spatialInertia(bodyMass, bodyCOM, bodyRotInertia);
bodyInertia = mcI(bodyMass, bodyCOM([1 3]), bodyRotInertia(2,2) );

abadLocation = [bodyLength, bodyWidth, 0] * 0.5;
hipLocation = [0, abadLinkLength, 0];
kneeLocation = [0, 0, -hipLinkLength];
footLocation = [0, 0, -kneeLinkLength];

hipSrbmLocation = [0.19 -0.1 0;0.19 0.1 0;...
    -0.19 -0.1 0;-0.19 0.1 0];

robot = struct(...
    'bodyMass',bodyMass,...
    'bodyLength', bodyLength,...
    'bodyWidth', bodyWidth,...
    'bodyHeight',bodyHeight,...
    'bodyInertia',bodyInertia,...
    'abadInertia',abadInertia,...
    'abadLocation',abadLocation,...
    'hipInertia',hipInertia,...
    'hipLocation',hipLocation,...
    'kneeInertia',kneeInertia,...
    'kneeLocation',kneeLocation,...
    'footLocation',footLocation,...
    'hipSrbmLocation',hipSrbmLocation,...
    'leg_rad',.02,...
    'abadGearRatio',6.0,...
    'hipGearRatio',6.0,...
    'kneeGearRatio',9.33,...
    'motorKT',0.05,...
    'motorR',0.173,...
    'motorTauMax',3.0,...
    'motorVelMax',240.0,...
    'batteryV',24);

%% Build the Model
model = buildQuadrupedModel2D(robot);
model.name    = 'MiniCheetah2D';
model.params  = robot; 
