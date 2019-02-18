%{
Authors: Josh Penner, John Matson and Aileen Maynard
Date: February 15th, 2019
Course: MECH 4630 - Statics and Dynamics
Set: T
%}

clear %clear any previous data to start

fail = -1; %testing for file open failure

%{
brief: find user designated file and load information into program
param:
return:
%}
prompt = 'Please enter file to open: ';
filename = input(prompt,'s'); %user inputs full file name
fid = fopen(filename,'r') ; %chosen text file stored into fid

while fid == fail %if no file is opened
    prompt2 = 'File not found, please re-enter file to open: ';
    filename = input(prompt2,'s');
    fid = fopen(filename,'r') ;
end

S = textscan(fid,'%s','Delimiter','\n');
S = S{1} ;


%{
Brief: Parse the file information into usuable arrays
Param:
Return:
%}
%Get the line numbers of each section
idxS1 = strfind(S, '[JOINT COORDINATES]');
idx1 = find(not(cellfun('isempty', idxS1)));
idxS2 = strfind(S, '[MEMBER JOINT CONNECTIVITY]');
idx2 = find(not(cellfun('isempty', idxS2)));
idxS3 = strfind(S, '[REACTIONS AT NODES]');
idx3 = find(not(cellfun('isempty', idxS3)));
idxS4 = strfind(S, '[EXTERNAL FORCES]');
idx4 = find(not(cellfun('isempty', idxS4)));
idxS5 = strfind(S, '[FORCE UNITS]');
idx5 = find(not(cellfun('isempty', idxS5)));

frewind(fid);

% Read data from joint section into structure
jointsection = textscan(fid, '%f %f %f ','headerLines',idx1,'delimiter','\t');
joint.index=cell2mat(jointsection(1));
joint.x=cell2mat(jointsection(2));
joint.y=cell2mat(jointsection(3));

frewind(fid);

% Read data from member section into structure
membersection = textscan(fid, '%d %d %d','headerLines',idx2,'delimiter','\t');
member.index = cell2mat(membersection(1));
member.jointA = cell2mat(membersection(2));
member.jointB = cell2mat(membersection(3));

frewind(fid);

% Read data from reaction section into structure
reactionsection = textscan(fid, '%d %d %c','headerLines',idx3,'delimiter','\t');
reaction.index = cell2mat(reactionsection(1));
reaction.node = cell2mat(reactionsection(2));
reaction.dir = cell2mat(reactionsection(3));

frewind(fid);

% Read data from external force section into structure
extfsection = textscan(fid, '%d %d %f %c','headerLines',idx4,'delimiter','\t');
extf.index = cell2mat(extfsection(1));
extf.joint = cell2mat(extfsection(2));
extf.force = cell2mat(extfsection(3));
extf.dir = cell2mat(extfsection(4));

frewind(fid);

% Read data from unit section into structure
units = cell2mat(textscan(fid, '%c','headerLines',idx5));

fclose(fid) ; % close file

%{
Brief: Error check: number of variables from file
Param:
Return:
%}
NJ = size(joint.index,1); % number of joints

NM = size(member.index,1); % number of members

NR = size(reaction.index,1); % number of reactions

NE = size(extf.index,1); % number of external forces

if 2*NJ ~= (NM+NR)
    message = '2*(#Joints) does not equal (#Members + #Reactions)';
    fprintf ('\nNJ = %d, NM = %d and NR = %d\n', NJ, NM, NR)
    error(message)
end

%{
Brief: Create the matrix for calculations
Param:
Return:
%}
E = zeros(NJ*2,1); % create E matrix for external forces

% fill E matrix with external force values
for i = 1:NE
    Eindex = extf.joint(i) * 2 + 1;
    if (extf.dir(i) == 'Y')
        Eindex = Eindex + 1 ;
    end
    E(Eindex)=extf.force(i);
end

M = zeros(NJ*2,NM+NR); % create M matrix for structure shape

% load M matrix with coefficients for member forces - joint A
joint_index = 0;
for mem_index = 1:NM 
    ydist = joint.y(member.jointB(mem_index) + 1) - joint.y(member.jointA(mem_index) + 1);
    xdist = joint.x(member.jointB(mem_index) + 1) - joint.x(member.jointA(mem_index) + 1);
    length = sqrt(ydist^2 + xdist ^2);
    Mcol = mem_index;
    joint_index = member.jointA(mem_index);
    Mrow = (joint_index * 2) + 1;
    M(Mrow,Mcol) = xdist / length;
    Mrow = Mrow + 1;
    M(Mrow,Mcol) = ydist / length;
    
    joint_index = member.jointB(mem_index);
    Mrow = (joint_index * 2) + 1;
    M(Mrow,Mcol) = -xdist / length;
    Mrow = Mrow + 1;
    M(Mrow,Mcol) = -ydist / length;
end

% load M matrix with coefficients for reaction forces
for reaction_index = 1:NR 
    Mcol = reaction_index + NM;
    Mrow = (reaction.node(reaction_index) * 2) + 1;
    if (reaction.dir(reaction_index) == 'Y')
        Mrow = Mrow + 1;
    end
    M(Mrow,Mcol) = 1;
end

result = inv(M) * (-E);

%{
Brief: Places matrix and results into user designated file
Param:
Return:
%}
text = '.txt';
prompt = 'Please enter file to write to: ';
filename = input(prompt,'s');
filename = strcat(filename,text);
trust = fopen(filename,'wt') ;
fprintf(trust,'M Matrix\n--------\n');

rowstring = "";

for row = 1:(NM+NR)
    for col = 1:(NM+NR)
           rowstring = strcat(rowstring, num2str(M(row,col), '%+.2f'), "\t");
     end
    fprintf(trust,rowstring);
    fprintf(trust,"\n");
    rowstring = "";     
end

fprintf(trust,'\n\nSolution\n--------\n');

counter = 0;
 
for memnum = 1:NM
    fprintf(trust,'Member %02.0f: F= %10.2f ', counter, abs(result(memnum)));
    fprintf(trust,units);
    if result(memnum)<0
        fprintf(trust, ' [C]');
    else
        fprintf(trust, ' [T]');
    end
    fprintf(trust, '\n');
    counter=counter + 1;
end

fprintf(trust, '\n');

for reactnum = 1:NR
    if reaction.dir(reactnum) == 'X'
        fprintf(trust, 'X');
    else
        fprintf(trust, 'Y');
    end 
    fprintf(trust,'-direction reaction on joint %2.0f = %+10.2f ', reaction.node(reactnum), result(NM+reactnum));
    fprintf(trust,units);
    fprintf(trust,'\n');
end

fclose(trust) ; % close file

