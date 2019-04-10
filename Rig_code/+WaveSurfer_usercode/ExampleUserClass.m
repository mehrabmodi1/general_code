classdef ExampleUserClass < ws.UserClass

    % This is a very simple user class.  It writes to the console when
    % things like a sweep start/end happen.
    
    % Information that you want to stick around between calls to the
    % functions below, and want to be settable/gettable from outside the
    % object.
    properties
        Greeting = 'Hello, there!'
        TimeAtStartOfLastRunAsString_ = ''  
          % TimeAtStartOfLastRunAsString_ should only be accessed from 
          % the methods below, but making it protected is a pain.
    end
    
    methods        
        function self = ExampleUserClass(parent)
            % creates the "user object"
            fprintf('%s  Instantiating an instance of ExampleUserClass.\n', ...
                    self.Greeting);
        end
        
        function delete(self)
            % Called when there are no more references to the object, just
            % prior to its memory being freed.
%             fprintf('%s  An instance of ExampleUserClass is being deleted.\n', ...
%                     self.Greeting);
        end
        
        % These methods are called in the WaveSurfer frontend process
        function startingRun(self,wsModel,eventName)
            % Called just before each set of sweeps (a.k.a. each
            % "run")
%             self.TimeAtStartOfLastRunAsString_ = datestr( clock() ) ;
%             fprintf('%s  About to start a run.  Current time: %s\n', ...
%                     self.Greeting,self.TimeAtStartOfLastRunAsString_);
%             
        end
        
        function completingRun(self,wsModel,eventName)
            % Called just after each set of sweeps (a.k.a. each
            % "run")
%             fprintf('%s  Completed a run.  Time at start of run: %s\n', ...
%                     self.Greeting,self.TimeAtStartOfLastRunAsString_);
        end
        
        function stoppingRun(self,wsModel,eventName)
            % Called if a sweep goes wrong
%             fprintf('%s  User stopped a run.  Time at start of run: %s\n', ...
%                     self.Greeting,self.TimeAtStartOfLastRunAsString_);
        end        
        
        function abortingRun(self,wsModel,eventName)
            % Called if a run goes wrong, after the call to
            % abortingSweep()
%             fprintf('%s  Oh noes!  A run aborted.  Time at start of run: %s\n', ...
%                     self.Greeting,self.TimeAtStartOfLastRunAsString_);
        end
        
        function startingSweep(self,wsModel,eventName)
            timestamp = now;
            % Called just before each sweep
%             fprintf('%s  About to start a sweep.  Time at start of run: %s\n', ...
%                     self.Greeting,self.TimeAtStartOfLastRunAsString_);
        end
        
        function completingSweep(self,wsModel,eventName)
            %checking if WaveSurfer is saving data
            if wsModel.IsLoggingEnabled == 1
                %identifying current aquisition directory
                curr_direc = curr_aq_direc();
                %checking if this is the first trial, loading in data if it
                %isn't
                if exist([curr_direc, '\axon_amp_metadata.mat']) == 2
                    metadata = load([curr_direc, '\axon_amp_metadata.mat']);
                    metadata = metadata.metadata;
                else
                    metadata = [];
                end

                a = ws.dabs.axon.MulticlampTelegraph('getAllElectrodeIDs');
                info = ws.dabs.axon.MulticlampTelegraph('getElectrodeState', a);
                info.timestamp = timestamp;           %adding on a timestamp

                curr_tr = size(metadata, 2) + 1;
                metadata(curr_tr).info = info;
                save([curr_direc, '\axon_amp_metadata.mat'], 'metadata');
            else
            end
            
            % Called after each sweep completes
%             fprintf('%s  Completed a sweep.  Time at start of run: %s\n', ...
%                     self.Greeting,self.TimeAtStartOfLastRunAsString_);
        end
        
        function stoppingSweep(self,wsModel,eventName)
            % Called if a sweep goes wrong
%             fprintf('%s  User stopped a sweep.  Time at start of run: %s\n', ...
%                     self.Greeting,self.TimeAtStartOfLastRunAsString_);
        end        
        
        function abortingSweep(self,wsModel,eventName)
            % Called if a sweep goes wrong
%             fprintf('%s  Oh noes!  A sweep aborted.  Time at start of run: %s\n', ...
%                     self.Greeting,self.TimeAtStartOfLastRunAsString_);
        end        
        
        function dataAvailable(self,wsModel,eventName)
            % Called each time a "chunk" of data (typically 100 ms worth) 
            % has been accumulated from the looper.
%             analogData = wsModel.getLatestAIData() ;
%             digitalData = wsModel.getLatestDIData() ; 
%             nScans = size(analogData,1);
%             fprintf('%s  Just read %d scans of data.\n',self.Greeting,nScans);                                    
        end
        
        % These methods are called in the looper process
        function samplesAcquired(self,looper,eventName,analogData,digitalData) 
            % Called each time a "chunk" of data (typically a few ms worth) 
            % is read from the DAQ board.
%             nScans = size(analogData,1);
%             fprintf('%s  Just acquired %d scans of data.\n',self.Greeting,nScans);                                    
        end
        
        % These methods are called in the refiller process
        function startingEpisode(self,refiller,eventName)
            % Called just before each episode
%             fprintf('%s  About to start an episode.\n',self.Greeting);
        end
        
        function completingEpisode(self,refiller,eventName)
            % Called after each episode completes
%             fprintf('%s  Completed an episode.\n',self.Greeting);
        end
        
        function stoppingEpisode(self,refiller,eventName)
            % Called if a episode goes wrong
%             fprintf('%s  User stopped an episode.\n',self.Greeting);
        end        
        
        function abortingEpisode(self,refiller,eventName)
            % Called if a episode goes wrong
%             fprintf('%s  Oh noes!  An episode aborted.\n',self.Greeting);
        end
    end  % methods
    
end  % classdef

