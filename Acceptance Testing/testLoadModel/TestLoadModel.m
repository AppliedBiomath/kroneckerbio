
classdef TestLoadModel < matlab.unittest.TestCase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        model
    end
    
    methods (TestMethodSetup)
        function loadModel(testCase)
            testCase.model = LoadModel('testModel.txt');
        end
    end
    
    methods (Test)
        function testLoadModelCompartments(testCase)
            % There should be two compartments "plasma" and "tissue"
            % check that they are 3dimentional and have the right names
            % and sizes
            c = testCase.model.Compartments;
            testCase.verifySize(c, [2,1]);
            
            testCase.verifyEqual(c(1).Size, 5);
            testCase.verifyEqual(c(1).Name, 'plasma');
            testCase.verifyEqual(c(1).Dimension, 3);
            
            testCase.verifyEqual(c(2).Size, 1);
            testCase.verifyEqual(c(2).Name, 'tissue');
            testCase.verifyEqual(c(2).Dimension, 3);
        end
        
        function testLoadModelStates(testCase)
            % There should be two States "a" and "b" and they should both
            % be in the plasma and have the correct initial conditions.
            s = testCase.model.States;
            testCase.verifySize(s, [6,1]);
            
            testCase.verifyEqual(s(1).Name, 'a');
            testCase.verifyEqual(s(1).Compartment, 'plasma');
            testCase.verifyEqual(s(1).InitialValue{2}, 5)
            
            testCase.verifyEqual(s(2).Name, 'b');
            testCase.verifyEqual(s(2).Compartment, 'plasma');
            testCase.verifyEqual(s(2).InitialValue{2}, 0);
            
            testCase.verifyEqual(s(3).Name, 'c');
            testCase.verifyEqual(s(3).Compartment, 'plasma');
            testCase.verifyEqual(s(3).InitialValue{1}, 'c0');
            testCase.verifyEqual(s(3).InitialValue{2}, 5);    
            
            testCase.verifyEqual(s(4).Name, 'a');
            testCase.verifyEqual(s(4).Compartment, 'tissue');
            testCase.verifyEqual(s(4).InitialValue{2}, 5)
            
            testCase.verifyEqual(s(5).Name, 'b');
            testCase.verifyEqual(s(5).Compartment, 'tissue');
            testCase.verifyEqual(s(5).InitialValue{2}, 0);
            
            testCase.verifyEqual(s(6).Name, 'c');
            testCase.verifyEqual(s(6).Compartment, 'tissue');
            testCase.verifyEqual(s(6).InitialValue{1}, 'c0');
            testCase.verifyEqual(s(6).InitialValue{2}, 1);       
        end
        
        function testLoadModelParameters(testCase)
            % There should be two rate constants "kfwd" and "kref" and they
            % should have the correct values.
            p = testCase.model.Parameters;
            
            testCase.verifySize(p, [2,1]);
            
        end
        
        function testLoadModelFlux(testCase)
           t  = 0;
           u0 = zeros(0, 1);
           f  = testCase.model.f;

           testCase.verifyEqual(f(t, [1; 1; 0; 1; 1; 0], u0), [0;0;0;0;0;0]);
           testCase.verifyEqual(f(t, [1; 0; 0; 0; 1; 0], u0), [-1;1;0;1;-1;0]);
        end
        
        
        function testLoadModelInitialCondion(testCase)
            m = testCase.model;
            con = Experiment(m, 1);

            % Run the simulations
            sim1 = Simulate(m, con);

            testIC = [5;0;50;5;0;10];
            testCase.verifyEqual(sim1.x(0), testIC);
        end
        
    end
    
end

