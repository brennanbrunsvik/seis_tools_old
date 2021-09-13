function [pz]=readsacpz_rdseed(varargin)
%READSACPZ_RDSEED    Reads in SAC polezero files generated by RDSEED
%
%    Usage:    pz=readsacpz_rdseed(file)
%              pz=readsacpz_rdseed(file1,...,fileN)
%              pz=readsacpz_rdseed(string,true)
%
%    Description:
%     PZ=READSACPZ_RDSEED(FILE) reads in an "annotated" SAC polezero file
%     generated by the program RDSEED (version 5.1+) and outputs a scalar
%     structure PZ containing all the information therein.  In particular,
%     this parses the information in the comment block above the pole,
%     zero, & constant info unlike READSACPZ which just discards it.
%     Multiple polezero entries in a single file are supported (as output
%     by the IRIS sacpz web service).  The output struct PZ will have
%     fields of size NPZx1 where NPZ is the number of polezero entries.
%
%     PZ=READSACPZ_RDSEED(FILE1,...,FILEN) reads in multiple SAC polezero
%     files and concatenates the resulting structures (if possible)
%     returning a single structure PZ.  If polezero files have differing
%     RDSEED comment block fields then concatenation will fail, a warning
%     is generated & the output is a cell array of structures that would be
%     returned as if READSACPZ_RDSEED was called on each file individually.
%
%     PZ=READSACPZ_RDSEED(STRING,TRUE) allows parsing a string as if it was
%     the contents of a SAC polezero file.  This is useful for parsing info
%     directly from the IRIS sacpz web service (see Examples section).
%
%    Notes:
%     - Each polezero MUST have a comment block to be parsed correctly.  If
%       one does not exist an error will be generated or if the polezero
%       comes after an annotated polezero in the same file the polezero
%       info will overwrite that previous annotated polezero's zeros, poles
%       and constant.
%     - Multiple ZEROS, POLES or CONSTANT lines after a comment block are
%       allowed but not encouraged as they overwrite the previous entry.
%       For instance, when there is two CONSTANT lines the first one will
%       be read and overwrote.
%     - An example of a RDSEED "annotated" PoleZero:
%        * **********************************
%        * NETWORK   (KNETWK): BK
%        * STATION    (KSTNM): ARC
%        * LOCATION   (KHOLE): 
%        * CHANNEL   (KCMPNM): BHE
%        * CREATED           : 2014-03-06T15:40:13
%        * START             : 1992-05-28T19:48:00
%        * END               : 2001-08-02T24:36:60
%        * DESCRIPTION       : Humboldt State University, Arcata, CA, USA
%        * LATITUDE          : 40.877720 
%        * LONGITUDE         : -124.077380
%        * ELEVATION         : 30.1   
%        * DEPTH             : 0.0  
%        * DIP               : 90.0 
%        * AZIMUTH           : 90.0 
%        * SAMPLE RATE       : 20.0
%        * INPUT UNIT        : M
%        * OUTPUT UNIT       : COUNTS
%        * INSTTYPE          : Streckeisen STS-2 VBB Tri-Axial Seismometer
%        * INSTGAIN          : 1.500000e+03 (M/S)
%        * COMMENT           : N/A
%        * SENSITIVITY       : 6.373550e+08 (M/S)
%        * A0                : 4.862460e+07
%        * **********************************
%        ZEROS	3
%        	+0.000000e+00	+0.000000e+00
%        	+0.000000e+00	+0.000000e+00
%        	+0.000000e+00	+0.000000e+00
%        POLES	5
%        	-3.702370e-02	+3.702440e-02
%        	-3.702370e-02	-3.702440e-02
%        	-1.187520e+02	+4.234880e+02
%        	-1.187520e+02	-4.234880e+02
%        	-2.513270e+02	+0.000000e+00
%        CONSTANT	+3.099113e+16
%
%    Examples:
%     % Read in all polezero info from station ANMO:
%     url=['http://service.iris.edu/irisws/sacpz/1/' ...
%          'query?net=IU&loc=*&cha=*&sta=ANMO'];
%     pz=readsacpz_rdseed(urlread(url),true);
%
%    See also: READSACPZ, WRITESACPZ_RDSEED, WRITESACPZ, REMOVESACPZ,
%              APPLYSACPZ, MAKESACPZDB, GENSACPZNAME, PARSE_SACPZ_FILENAME,
%              GETSACPZ, ISSACPZ_RDSEED, FIX_OLD_SACPZ, SSIDX, SSCAT

%     Version History:
%        Feb. 25, 2014 - initial version
%        Mar.  5, 2014 - minor doc update, input/output handles multi-word
%        Mar.  6, 2014 - doc update, bugfix: allow CONSTANT to be NaN
%        Apr. 28, 2014 - bugfix: preallocation fails no more if no ZEROS
%        May  28, 2014 - bugfix: unconcatenatable output numeric fields are
%                        now converted from strings, trim excess allocation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  28, 2014 at 15:05 GMT

% todo:

% string input or files?
% addpath('/Users/zeilon/Documents/MATLAB/seizmo/time')
global SEIZMO
if(nargin==2 && islogical(varargin{2}) && isscalar(varargin{2}) ...
        && varargin{2} && isstring_seizmo(varargin{1}))
    nfiles=1;
    string=true;
else
    % files
    string=false;
    
    % set filterspec appropriately
    SEIZMO.ONEFILELIST.FILTERSPEC=...
        {'SAC_PZ*;sac_pz*;SACPZ*;sacpz*' ...
        'SACPZ Files (SAC_PZ*,sac_pz*,sacpz*,SACPZ*)'};
    
    % compile file lists
    varargin=onefilelist(varargin{:});
    nfiles=numel(varargin);
    
    % error if no files
    if(nfiles<1)
        error('seizmo:readsacpz_rdseed:noFilesToRead','No files to read!');
    end
end


% detail message
verbose=false;
if(verbose)
    disp('Reading In SAC Polezero File(s)');
    print_time_left(0,nfiles);
end

% loop over each file
pz=cell(nfiles,1);
for f=1:nfiles
    % string input
    if(string)
        a=getwords(varargin{1},sprintf('\n'));
        clear varargin;
        varargin(f).path=['.' filesep];
        varargin(f).name='SAC_PZs_FILENAME_UNKNOWN';
    else
        % read in text from file
        file=[varargin(f).path varargin(f).name];
        a=getwords(readtxt(file),sprintf('\n'));
    end

    % error on empty file
    if(isempty(a))
        error('seizmo:readsacpz_rdseed:emptySACPZ',...
            'SAC PoleZero File: %s\nFile is empty!',file);
    end
    
    % estimate number of polezero entries in this file
    npz_guess=max(1,sum(strncmp('ZEROS',a,5)));
    pre(1:npz_guess,1)={''}; % preallocation
    
    % keep processing until through all lines
    npz=0; line=1; nlines=numel(a);
    while(line<=nlines)
        % skip line if blank
        if(isempty(a{line}))
            line=line+1;
            continue;
        end
        
        % process line
        words=getwords(a{line});
        
        % skip line if blank or comment unless it is the beginning
        % of a rdseed comment header that we want to parse
        if(isempty(words) || strcmp(words{1}(1),'*'))
            if(~isempty(words) && strcmp(words{1},'*') ...
                    && numel(words{2})>=10 ...
                    && strcmp(words{2}(1:10),'**********'))
                % rdseedv5+ comment block header
                npz=npz+1;
                if(npz==1)
                    pz{f}.path=pre;
                    pz{f}.path=pre;
                end
                pz{f}.path{npz,1}=varargin(f).path;
                pz{f}.name{npz,1}=varargin(f).name;
                
                % parse lines until end of comment block
                line=line+1;
                while(line<=nlines)
                    % process line
                    words=getwords(a{line});
                    if(isempty(words))
                        % empty line in comment block...user modification?
                        line=line+1;
                    elseif(strcmp(words{1},'*') && numel(words{2})>=10 ...
                            && strcmp(words{2}(1:10),'**********'))
                        % end of rdseedv5+ comment block
                        break;
                    elseif(strcmp(words{1},'*'))
                        % comment block info to parse
                        switch lower(words{2})
                            case 'network'
                                % use sac field
                                if(npz==1); pz{f}.knetwk=pre; end
                                pz{f}.knetwk(npz,1)=words(4);
                            case 'station'
                                % use sac field
                                if(npz==1); pz{f}.kstnm=pre; end
                                pz{f}.kstnm(npz,1)=words(4);
                            case 'location'
                                % use sac field
                                if(npz==1); pz{f}.khole=pre; end
                                if(numel(words)==4)
                                    pz{f}.khole(npz,1)=words(4);
                                else
                                    % empty location code
                                    pz{f}.khole{npz,1}='';
                                end
                            case 'channel'
                                % use sac field
                                if(npz==1); pz{f}.kcmpnm=pre; end
                                pz{f}.kcmpnm(npz,1)=words(4);
                            case 'start'
                                % convert later together
                                if(npz==1); pz{f}.b=pre; end
                                pz{f}.b(npz,1)=words(4);
                            case 'end'
                                % convert later together
                                if(npz==1); pz{f}.e=pre; end
                                pz{f}.e(npz,1)=words(4);
                            case 'latitude'
                                % use sac field
                                if(npz==1); pz{f}.stla=pre; end
                                pz{f}.stla(npz,1)=words(4);
                            case 'longitude'
                                % use sac field
                                if(npz==1); pz{f}.stlo=pre; end
                                pz{f}.stlo(npz,1)=words(4);
                            case 'elevation'
                                % use sac field
                                if(npz==1); pz{f}.stel=pre; end
                                pz{f}.stel(npz,1)=words(4);
                            case 'depth'
                                % use sac field
                                if(npz==1); pz{f}.stdp=pre; end
                                pz{f}.stdp(npz,1)=words(4);
                            case 'dip'
                                % use sac field
                                if(npz==1); pz{f}.cmpinc=pre; end
                                pz{f}.cmpinc(npz,1)=words(4);
                            case 'azimuth'
                                % use sac field
                                if(npz==1); pz{f}.cmpaz=pre; end
                                pz{f}.cmpaz(npz,1)=words(4);
                            case 'sample' % RATE
                                if(npz==1); pz{f}.sr=pre; end
                                pz{f}.sr(npz,1)=words(5);
                            case 'input' % UNIT
                                if(npz==1); pz{f}.input=pre; end
                                if(numel(words)>=5)
                                    pz{f}.input{npz,1}=...
                                        joinwords(words(5:end));
                                else
                                    % maybe input==output?
                                    pz{f}.input{npz,1}='';
                                end
                            case 'output' % UNIT
                                if(npz==1); pz{f}.output=pre; end
                                if(numel(words)>=5)
                                    pz{f}.output{npz,1}=...
                                        joinwords(words(5:end));
                                else
                                    % maybe input==output?
                                    pz{f}.output{npz,1}='';
                                end
                            case 'instgain'
                                if(npz==1); pz{f}.instgain=pre; end
                                if(npz==1); pz{f}.instgainunits=pre; end
                                if(numel(words)>3)
                                    pz{f}.instgain(npz,1)=words(4);
                                else
                                    pz{f}.instgain{npz,1}='';
                                end
                                if(numel(words)>4)
                                    pz{f}.instgainunits(npz,1)=words(5);
                                else
                                    pz{f}.instgainunits{npz,1}='';
                                end
                            case 'sensitivity'
                                if(npz==1); pz{f}.sensitivity=pre; end
                                if(numel(words)>3)
                                    pz{f}.sensitivity(npz,1)=words(4);
                                else
                                    pz{f}.sensitivity{npz,1}='';
                                end
                                if(numel(words)>4)
                                    pz{f}.sensitivityunits(npz,1)=words(5);
                                else
                                    pz{f}.sensitivityunits{npz,1}='';
                                end
                            case 'a0'
                                if(npz==1); pz{f}.a0=pre; end
                                pz{f}.a0(npz,1)=words(4);
                            otherwise
                                if(npz==1)
                                    pz{f}.(lower(words{2}))=pre;
                                end
                                pz{f}.(lower(words{2})){npz,1}=...
                                    joinwords(words(4:end));
                        end
                        line=line+1;
                    else
                        % no idea - this shouldn't happen unless
                        %           the file was modified outside rdseed
                        if(strcmp(words{1}(1),'*'))
                            % user comment
                            line=line+1;
                        else
                            % user deleted part of rdseed comment block?
                            % Regardless, this is a noncomment.  So
                            % compensate for the line incrementing later
                            % on and break out of the comment block loop.
                            warning(...
                                'seizmo:readsacpz_rdseed:corruptSACPZ',...
                                ['SAC PoleZero File: %s\n' ...
                                'Found corrupted comment block!'],file);
                            line=line-1;
                            break;
                        end
                    end
                end
                
                % default z/p/k (gives unity at all frequencies
                if(npz==1)
                    pz{f}.z=cell(npz_guess,1);
                    pz{f}.p=cell(npz_guess,1);
                    pz{f}.k=ones(npz_guess,1);
                end
            end
            line=line+1;
            continue;
        end
        
        % check that a good comment block was found
        if(npz==0 || isempty(pz{f}.kstnm{npz,1}))
            error('seizmo:readsacpz_rdseed:notRdseedSACPZ',...
                'No RDSEED SACPZ comment block found! Use READSACPZ!');
        end
        
        % check number of words
        if(numel(words)~=2)
            error('seizmo:readsacpz_rdseed:notSACPZ',...
                ['SAC PoleZero File: %s\nLine %d: %s\n'...
                'Does Not Conform To SAC PoleZero Format!'],...
                file,line,a{line});
        end
        
        % separate words
        [field,value]=deal(words{:});
        value=str2double(value);
        
        % check value
        if(isnan(value) && ~strcmpi('CONSTANT',field))
            error('seizmo:readsacpz_rdseed:notSACPZ',...
                ['SAC PoleZero File: %s\nLine %d: %s\n'...
                'Does Not Conform To SAC PoleZero Format!'],...
                file,line,a{line});
        end
        
        % act based on line
        switch lower(field)
            case {'zeros' 'poles'}
                % check value is fixed
                if(value~=fix(value))
                    error('seizmo:readsacpz_rdseed:notSACPZ',...
                        ['SAC PoleZero File: %s\nLine %d: %s\n'...
                        'Does Not Conform To SAC PoleZero Format!'],...
                        file,line,a{line});
                end
                
                % preallocate poles/zeros
                switch lower(field)
                    case 'zeros'
                        pz{f}.z{npz,1}=zeros(value,1);
                    case 'poles'
                        pz{f}.p{npz,1}=zeros(value,1);
                end
                
                % read in poles/zeros (not at origin)
                line=line+1; nv=1;
                while(line<=nlines && nv<=value)
                    % skip line if blank
                    if(isempty(a{line}))
                        line=line+1;
                        continue;
                    end
                    
                    % process line
                    words=getwords(a{line});
                    
                    % skip line if blank or comment
                    if(isempty(words) || strcmp(words{1}(1),'*'))
                        line=line+1;
                        continue;
                    end
                    
                    % check number of words
                    if(numel(words)~=2)
                        error('seizmo:readsacpz_rdseed:notSACPZ',...
                            ['SAC PoleZero File: %s\nLine %d: %s\n'...
                            'Does Not Conform To SAC PoleZero Format!'],...
                            file,line,a{line});
                    end
                    
                    % separate words
                    [value1,value2]=deal(words{:});
                    value1=str2double(value1);
                    value2=str2double(value2);
                    
                    % check values
                    if(isnan(value1) || isnan(value2)); break; end
                    
                    % assign values
                    switch lower(field)
                        case 'zeros'
                            if(value2) % complex
                                pz{f}.z{npz,1}(nv)=value1+value2*1i;
                            else % real
                                pz{f}.z{npz,1}(nv)=value1;
                            end
                        case 'poles'
                            if(value2) % complex
                                pz{f}.p{npz,1}(nv)=value1+value2*1i;
                            else % real
                                pz{f}.p{npz,1}(nv)=value1;
                            end
                    end
                    
                    % increment
                    nv=nv+1;
                    line=line+1;
                end
            case 'constant'
                pz{f}.k(npz,1)=value;
                line=line+1;
            otherwise
                error('seizmo:readsacpz_rdseed:notSACPZ',...
                    ['SAC PoleZero File: %s\nLine %d: %s\n'...
                    'Does Not Conform To SAC PoleZero Format!'],...
                    file,line,a{line});
        end
    end
    
    % trim off excess preallocation
    pz{f}=ssidx(pz{f},1:npz);
    
    % detail message
    if(verbose); print_time_left(f,nfiles); end
end

% attempt to concatenate
try
    pz=sscat(pz{:});
catch
    warning('seizmo:readsacpz_rdseed:multiheader',...
        ['The comment block table has changed between files.  Output ' ...
        'will be a cell-array of SAC polezero structs!']);
end

% detail message
if(verbose)
    disp('Converting SAC Polezero Strings to Numeric Values');
end

% now fix the rdseedv5+ header fields
if(~iscell(pz))
    pz=fix_rdseed_hdr(pz);
else
    for i=1:numel(pz)
        pz{i}=fix_rdseed_hdr(pz{i});
    end
end

end


function [pz]=fix_rdseed_hdr(pz)
% converts comment block character strings to a better format
tmp=datevec(pz.created,'yyyy-mm-ddTHH:MM:SS');
pz.created=[cal2doy(tmp(:,1:3)) tmp(:,4:end)];
tmp=datevec(pz.b,'yyyy-mm-ddTHH:MM:SS');
pz.b=[cal2doy(tmp(:,1:3)) tmp(:,4:end)];
tmp=datevec(pz.e,'yyyy-mm-ddTHH:MM:SS');
pz.e=[cal2doy(tmp(:,1:3)) tmp(:,4:end)];
pz.stla=str2double(pz.stla);
pz.stlo=str2double(pz.stlo);
pz.stel=str2double(pz.stel);
pz.stdp=str2double(pz.stdp);
pz.cmpinc=str2double(pz.cmpinc);
pz.cmpaz=str2double(pz.cmpaz);
pz.sr=str2double(pz.sr);
pz.instgain=str2double(pz.instgain);
pz.sensitivity=str2double(pz.sensitivity);
pz.a0=str2double(pz.a0);
end
