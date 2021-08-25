% load binary gdf file in structarray
%
% "gdf" is a cell array with data
%     with two fields:
%           .data => structure with data in single group
%           .params => structure with params of that group
% "info" is structured array with main header information

% Note: 
% 'fread' should be used with 'uint8=>char' instead of '*char', to prevent
% issues with utf-8 characters

function [gdf, info]= loadgdf(filename,arrays)

% Check number of arguments
if( nargin==1 ) 
    arrays = []; 
end

% Constants
GDFNAMELEN = 16;                 % Length of the ascii-names
GDFID  = 94325877;               % ID for GDF      

% Data types 
t_ascii  = hex2dec('0001');      % ASCII character
t_s32    = hex2dec('0002');      % Signed long           
t_dbl    = hex2dec('0003');      % Double    
t_undef  = hex2dec('0000');      % Data type not defined
t_null	 = hex2dec('0010');      % No data
t_u8	 = hex2dec('0020');      % Unsigned char
t_s8	 = hex2dec('0030');      % Signed char
t_u16	 = hex2dec('0040');      % Unsigned short
t_s16	 = hex2dec('0050');      % Signed short
t_u32	 = hex2dec('0060');      % Unsigned long
t_u64	 = hex2dec('0070');      % Unsigned 64bit int
t_s64	 = hex2dec('0080');      % Signed 64bit int
t_flt	 = hex2dec('0090');      % Float

% Block types 
t_dir    = hex2dec('0100');      % Directory entry start
t_edir   = hex2dec('0200');      % Directory entry end 
t_sval   = hex2dec('0400');      % Single valued         
t_arr    = hex2dec('0800');      % Array                 

% Open file
fid = fopen( filename );
if( fid < 0 ) 
    error('Error: file not found.'); 
end

% Read GDF main header
info = struct;

ID = fread(fid, 1, '*uint32');
if( ID ~= GDFID ) 
    error('Error: this is not a gdf file!') 
end

info.time_created = fread( fid, 1, '*uint32' );

creator = fread( fid, GDFNAMELEN, 'uint8=>char' )';
info.creator = creator( 1:find( creator==0, 1 )-1 );                % get string part upto zero-character

dest = fread( fid, GDFNAMELEN, 'uint8=>char' )';
info.destination = dest( 1:find( dest==0, 1 )-1 );                  % get string part upto zero-character

major = fread( fid, 1, '*uint8' );
minor = fread( fid, 1, '*uint8' );
info.gdf_version = strcat( num2str( major ), '.', num2str( minor ) );

major = fread( fid, 1, '*uint8' );
minor = fread( fid, 1, '*uint8' );
info.creator_version = strcat( num2str( major ), '.', num2str( minor ) );

major = fread( fid, 1, '*uint8' );
minor = fread( fid, 1, '*uint8' );
info.destination_version = strcat( num2str( major ), '.', num2str( minor ) );

fread( fid, 2, '*uint8' ); % read to next block

% Read GDF data blocks
gdf   = cell(1,1);
param = struct;
dat   = struct;
groupcount = 1;
level = 1;
lastarr = false;
while( ~feof(fid) )    
    % Read GDF block header
    name    = fread(fid, GDFNAMELEN, 'uint8=>char')';
    type    = fread(fid, 1, '*uint32');                             % block type
    size    = fread(fid, 1, '*uint32');                             % byte count of data in block
    
    % Get name
    name    = name( 1:find( name==0,1 )-1 );                        % get string part upto zero-character
    name    = genvarname(name);                                     % get rid of not allowed characters etc
    
    % Get block type
    dir  = ( bitand(type,t_dir)  >0 );
    edir = ( bitand(type,t_edir) >0 );
    sval = ( bitand(type,t_sval) >0 );
    arr  = ( bitand(type,t_arr)  >0 );
    
    % Get data type
    dattype = bitand(type, 255);
    
    % Check if array block is finished
    if( lastarr && ( isempty(arr) || ~arr ) )
          % Save data group
           gdf{groupcount,1}.data = dat;
           gdf{groupcount,1}.param = param;
           dat = struct;                                            % clear data arrays
           groupcount = groupcount + 1;                             % next data group         
    end
        
    % New folder
    if( dir ) 
        level = level + 1; 
    end
    
    % End folder 
    if( edir )                                        
        fields = fieldnames(param);
        param = rmfield(param, fields(end));                        % clear last added param                      
        level = level - 1;        
    end
        
    % Read single value
    if( sval )       
        switch( dattype )
            case t_dbl
                value = fread(fid, 1, 'double');                    %read param
                param.(name) = value;                               %save param
            case t_null
                % no data present
            case t_ascii
                value  = fread(fid, double(size), 'uint8=>char')';  %read param
                param.(name) = value;                               %save param
            case t_s32
                value  = fread(fid, 1, 'int32')';                   %read param
                param.(name) = value;                               %save param                
            otherwise
                disp('unknown datatype of value!!!');
                disp([ 'name=' name ]);
                disp([ 'type=' num2str(type,'%x') ]);
                disp([ 'size=' num2str(double(size),'%x') '\n' ]);                
                value = fread(fid, double(size), 'uint8=>char');    %skip data
        end
    end    
    
    % Read data array
    if( arr )          
        % check if array should be saved    
        if ( isempty(arrays)) || ~isempty( strcmp(name,arrays) )  
            savearray = true;
        else
            savearray = false;
        end
            
        switch( dattype )
            case t_dbl
                if( mod(size,8) ~= 0 )
                     error('error, wrong size for double array'); 
                end
                N = double( size/8 );                              %number of array elements
                value = fread(fid, N, 'double');                   %read array
                if( savearray )
                    dat.(name) = value;                            %save array
                end        
            otherwise
                disp('unknown datatype of array!!!');
                disp([ 'name=' name ]);
                disp([ 'type=' num2str(type,'%x') ]);
                disp([ 'size=' num2str(double(size),'%x') '\n' ]);
                value = fread(fid, double(size), 'uint8=>char');   %skip data
        end
        
    end
      
   lastarr = arr;   % save if last data block was an array (or not).
end

% Close file
fclose(fid);


