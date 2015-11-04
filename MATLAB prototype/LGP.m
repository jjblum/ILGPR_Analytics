classdef LGP < handle
    
    %{
        Represents one gaussian process.
        Uses GPML library
    %}
    
    properties
        X        % data locations matrix
        Z        % data value vector
        W        % length scale matrix
        hyp      % vector of log(hyperparameters)
        u        % center location
        L        % lower triangular (Cholesky) decomp of K
        LT       % transpose of lower triangular decomp of K (saved so you don't need to recompute)
        R        % permutation matrix
        K        % covariance matrix
        N        % number of training points
        NMAX     % maximum number of training points
        NSTART   % number of points needed before computation actually begins
        started  % boolean, false while N < NSTART
        sp       % activation
        mc       % mixture coefficient (pi in the paper)
        data     % vector of Datum objects
        f        % most recent likelihood result
        ID       % the index of the LGP itself (order of creation)
    end
    
    methods
        function obj = LGP(ID,Datum)
            
%             keyboard
            obj.ID = ID;
            obj.X = Datum.getX();
            obj.Z = Datum.getZ();
            obj.hyp = [log(5) log(5) log(1)]; % log(length scale 1), log(length scale 2), log(process variability)
            obj.W = diag((1./exp(obj.hyp(1:2))).^2);
            obj.u = Datum.getX();
            obj.R = [];
%             obj.K = covSEard(obj.hyp,obj.X'); % K = covSEard(hyp, x, z, i) - note that each datum x is a ROW, not a column
%             obj.L = chol(obj.K,'lower');
            %                 obj.K = [];
            %                 obj.L = [];
            obj.LT = [];
            obj.N = 1;
            obj.NMAX = 50;
            obj.NSTART = 5;
            obj.started = 0;
            obj.sp = 1;
            obj.data = cell(1,1);
            obj.data{1} = Datum;
            obj.mc = 0;
            obj.f = 0;
        end
        
        function newDatum(self,Datum)
            self.N = self.N + 1;
            self.data{self.N} = Datum;
            updateX(self,Datum.getX());
            updateZ(self,Datum.getZ());
            choleskyUpdate(self);
        end
        
        function updateX(self, x)
            self.X = horzcat(self.X,x);
            
        end
        
        function updateZ(self, z)
            self.Z = vertcat(self.Z,z);
        end
        
        function getR(self, m)
            % m = index of value being replaced in covariance matrix
        end
        
        function updateCenter(self, x, posterior)
%             keyboard
            self.u = self.u + posterior/self.sp*(x - self.u);
        end
        
        function choleskyUpdate(self)
            
%             keyboard
            
            if self.N < self.NSTART
                return;
            elseif self.N == self.NSTART
                firstCholesky(self);
            elseif self.N > self.NMAX
                % the incremental cholesky update with removal
            else
                % the incremental cholesky update without removal
            end
        end
        
        function incrementalCholeskyUpdate(self)
            
        end
        
        function firstCholesky(self)
%             keyboard
            firstK(self);
            self.L = chol(self.K,'lower');
            self.LT = self.L';
        end
        
        function firstK(self)
            keyboard
            self.K = covSEard(self.hyp,self.X');
        end
        
        function updateSP(self, p)
            self.sp = self.sp + p;
        end
        
        function sp = getSP(self)
            sp = self.sp;
        end
        
        function setMC(self,mc)
            self.mc = mc;
        end
        
        function mc = getMC(self)
            mc = self.mc;
        end
        
        function f = getF(self)
            f = self.f;
        end
        
        function updateF(self,x)
%             keyboard
            % ratio = norm(x-self.u)*sqrt(self.W(1,1))
            % rough_f = exp(-1/2*ratio)
            self.f = exp(-1/2*(x-self.u)'*self.W*(x-self.u)); 
            % approximately exp(-1/2*[ratio of squared distance to squared length scale]). e.g. if ratio is 5, f ~ 0.1
            % So if we want a new LGP when more than 2*length scale away, need to set cutoff to about 0.4
        end
    end
    
end