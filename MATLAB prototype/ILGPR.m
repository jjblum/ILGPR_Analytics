classdef ILGPR < handle
    
    %{
        
    %}
    
    properties
        LGPs % vector of LGP objects
        nLGPs % number of LGPs
        predictionX % grid of prediction locations
        predictionZ % avg. at prediction X
        predictionS % std dev. at prediction X
        newLGPCutoff % likelihood below which a new LGP will be added
        M % the maximum number of LGP models used in the prediction
        fSum % sum of the weights
    end
    
    methods
        function obj = ILGPR(predictionX, predictionZ, predictionS)
            obj.nLGPs = 0;
            obj.LGPs = cell(1,1);
            obj.predictionX = predictionX;
            obj.predictionZ = predictionZ;
            obj.predictionS = predictionS;
            obj.newLGPCutoff = 0.1;
            obj.M = 10;
            obj.fSum = 0;
        end       
        
        function newDatum(self,Datum)
            if self.nLGPs > 0 % there is at least one LGP
                likelihoods = zeros(self.nLGPs,1);
                
                % calculate likelihood for all LGPs
                for i = 1:self.nLGPs
                    self.LGPs{i}.updateF(Datum.getX());
                    likelihoods(i) = self.LGPs{i}.getF();
                end                
                
                [maxLikelihood,bestLGPIdx] = max(likelihoods);
                                
                if maxLikelihood > self.newLGPCutoff                                        
                    
                    self.LGPs{bestLGPIdx}.newDatum(Datum);
                    
                else % the max likelihood is too small, need new LGP 
                    addLGP(self,Datum);
                end
                
            else % the first LGP
                addLGP(self,Datum);
            end            
        end
        
        function addLGP(self,Datum)
            self.nLGPs = self.nLGPs + 1;
            newLGP = LGP(self.nLGPs,Datum);
            self.LGPs{self.nLGPs} = newLGP;
        end
        
        function [weightedZ,weightedS] = predict(self,x)
            
            numStarted = 0;
            validLGPIdx = [];
            for j=1:self.nLGPs
                if self.LGPs{j}.isStarted()
                    numStarted = numStarted + 1;
                    validLGPIdx = [validLGPIdx j];
                end
            end

            if numStarted > self.M
                numStarted = self.M;
            end
            if ~numStarted
                disp('WARNING: No started LGPs available for prediction');
                return;
            end
            
            self.fSum = 0;
            for j = validLGPIdx              
                self.LGPs{j}.updateF(x);
                self.fSum = self.fSum + self.LGPs{j}.getF();
            end                
            
            weights = zeros(numStarted,1);
            for j = 1:numStarted
                weights(j) = self.LGPs{validLGPIdx(j)}.getF()/self.fSum;
            end
            [sortedWeights, sorted_idx] = sort(weights,'descend');
%             sorted_idx = sorted_idx(1:numStarted);            
            


            weightedZ = 0;
            weightedS = 0;
            for j = sorted_idx'
%                 if self.LGPs{j}.isStarted() % no need to calculate unless it is started (because weight would be zero anyway)
                    [zHat,sHat] = self.LGPs{validLGPIdx(j)}.predict(x);
                    weightedZ = weightedZ + zHat*weights(j);
                    weightedS = weightedS + sHat*weights(j);
%                 end
            end
        end        
             
     
    end
    
end