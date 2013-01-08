%slra_obj Example MATLAB class wrapper to an underlying C++ class
classdef slra_obj < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance 
        function this = slra_obj(varargin)
            this.objectHandle = slra_mex_obj('new', varargin{:});
        end
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            slra_mex_obj('delete', this.objectHandle);
        end
        function varargout = optimizeGsl(this, varargin)
            [varargout{1:nargout}] = slra_mex_obj('optimize', this.objectHandle, varargin{:});
        end
        function varargout = getRini(this, varargin)
            [varargout{1:nargout}] = slra_mex_obj('getRini', this.objectHandle, varargin{:});
        end
        function varargout = getM(this, varargin)
            [varargout{1:nargout}] = slra_mex_obj('getM', this.objectHandle, varargin{:});
        end
        function varargout = func(this, varargin)
            [varargout{1:nargout}] = slra_mex_obj('func', this.objectHandle, varargin{:});
        end
        function varargout = grad(this, varargin)
            [varargout{1:nargout}] = slra_mex_obj('grad', this.objectHandle, varargin{:});
        end
        function varargout = getPh(this, varargin)
            [varargout{1:nargout}] = slra_mex_obj('getPh', this.objectHandle, varargin{:});
        end
    end
end
