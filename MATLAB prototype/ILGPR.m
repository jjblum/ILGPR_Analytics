classdef ILGPR < handle
    
    %{
        
    %}
    
    properties
        LGPs % vector of LGP objects
        nLGPs % number of LGPs
        predictionX % grid of prediction locations
        predictionZ % avg. at prediction X
        predictionS % std dev. at prediction X
        spSum % sum of LGP's activations
        mcfSum % sum of mixture components
        newLGPCutoff % likelihood below which a new LGP will be added
    end
    
    methods
        function obj = ILGPR(predictionX, predictionZ, predictionS)            
            obj.nLGPs = 0;
            obj.LGPs = cell(1,1);
            obj.predictionX = predictionX;
            obj.predictionZ = predictionZ;
            obj.predictionS = predictionS;
            obj.spSum = 0;
            obj.newLGPCutoff = 0.1;
        end
        
        function newDatum(self,Datum)
            
        end
        
        function findBestLGP(self,Datum)
            if self.nLGPs > 0
                likelihoods = zeros(self.nLGPs,1);
                for i = 1:self.nLGPs
                    self.LGPs{i}.updateF(Datum.getX());
                    likelihoods(i) = self.LGPs{i}.getF();
                end
                
                if min(likelihoods) > self.newLGPCutoff
                    updateMCFSum();
                    for i = 1:self.nLGPs
                        posterior = self.LGPs{i}.getMC()*self.LGPs{i}.getF()/self.mcfSum;
                        self.LGPs{i}.updateSP(posterior);
                    end
                else
                    addLGP(Datum);
                end
                
            else
                addLGP(Datum);
            end
            updateMCs();
        end
        
        function addLGP(self,Datum)
            self.nLGPs = self.nLGPs + 1;
            
            
        end        
        
        function prediction(self)
            
        end
        
        function updateMCs(self)
            updateSPSum();
            for i = 1:self.nLGPs
                self.LGPs{i}.setMC(self.LGPs{i}.getSP()/self.spSum);
            end               
        end
        
        function updateSPSum(self)
            self.spSum = 0;
            for i = 1:size(self.LGPs,1)
                self.spSum = self.spSum + self.LGPs{i}.getSP();
            end            
        end
        
        function updateMCFSum(self, x)
            self.mcfSum = 0;
            for i = 1:size(self.LGPs,1)
                self.mcfSum = self.LGPs{i}.getMC()*self.LGPs{i}.getF();
            end
        end
        

        
    end
    
end