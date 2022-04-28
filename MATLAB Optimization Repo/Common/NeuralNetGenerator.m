classdef NeuralNetGenerator
    % this class loads a ONNX file saved from pytorch and creates a matlab
    % neural network.
    properties
        matlab_net
    end
    methods
        function obj = NeuralNetGenerator(file)
            obj.matlab_net = importONNXNetwork(file,outputlayertype='regression');
        end
        function casadi_net = to_casadi(obj,use_jit)
            % converts neural net to casadi function. Warning: only works
            % for feedforward, fully connected.
            % optional argument : use_jit: in time compilation (true or
            % {false})
            
            if nargin > 1
                use_jit = use_jit;
            else
                use_jit=false;
            end
            
            if use_jit
                p_opts = append_jit_options(struct(),'-O3');
            else 
                p_opts =struct();
            end
            
            net = obj.matlab_net;
            in = casadi.MX.sym('x',net.Layers(1).InputSize,1);
            out = obj.eval_net(in);
            casadi_net = casadi.Function('net',{in},{out},p_opts);
            %%
            in_val = rand(net.Layers(1).InputSize,1);
            err = full(casadi_net(in_val)) - net.predict(in_val);
            
            tol=1e-6;
            if abs(err) > tol
                error(['Evualation test failed: error greater than ',num2str(tol)])
            else
                disp('Neural net imported successfully')
            end
        end
        
        function out = eval_net(obj,in)
            % Evaluates the layers of the matlab neural net. Only works for
            % fully connected NN. 
            
            % inputs : in : inputs to the neural network
            % outputs :out : outputs to the neural network
            
            net= obj.matlab_net;
            if length(in) ~= net.Layers(1).InputSize
                error('Wrong input size')
            end
            
            out=in;
            
            
            for i = 2:length(net.Layers)-1
                if contains(net.Layers(i).Name,'MatMul_')
                    out = double(net.Layers(i).Weights)*out + double(net.Layers(i).Bias);
                elseif contains(net.Layers(i).Name,'Add_')
                    out = out + double(net.Layers(i).Offset);
                elseif contains(net.Layers(i).Name,'Tanh_')
                    out=tanh(out);
                else
                    error('layer type not implemented yet')
                end
            end
        end
        
    end
end