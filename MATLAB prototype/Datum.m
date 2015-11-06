classdef Datum < handle
    
    %{
        
    %}
    
    properties
        x % position
        z % value
        i % unique ID
    end
    
    methods
        function obj = Datum(x,z,i)
            obj.x = x;
            obj.z = z;
            obj.i = i;
        end
        
        function x = getX(self)
            x = self.x;
        end
        
        function z = getZ(self)
            z = self.z;
        end
    end
    
end