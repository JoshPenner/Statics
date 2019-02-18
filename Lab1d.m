%size of data error check

clear %clear any previous data
fId = -1;

prompt = 'Please enter file to open: ';
fileName = input(prompt,'s');

fId = fopen(fileName,'r') ;

%{
while fid == -1
    prompt2 = 'File not found, please re-enter file to open: ';
    filename = input(prompt2,'s');
    fid = fopen(filename,'r') ;
end
%}

S = textscan(fId,'%s','Delimiter','\n');
S = S{1} ;

%%Get the line numbers of each section
idx1 = find(contains(S,'[JOINT COORDINATES]'));
idx2 = find(contains(S,'[MEMBER JOINT CONNECTIVITY]'));
idx3 = find(contains(S,'[REACTIONS AT NODES]'));
idx4 = find(contains(S,'[EXTERNAL FORCES]'));
idx5 = find(contains(S,'[FORCE UNITS]'));

frewind(fId);

% Read data from joint section into structure
jointSection = textscan(fId, '%f %f %f ','headerLines',idx1,'delimiter','\t');
joint.index=cell2mat(jointSection(1));
joint.x=cell2mat(jointSection(2));
joint.y=cell2mat(jointSection(3));

frewind(fId);

% Read data from member section into structure
memberSection = textscan(fId, '%d %d %d','headerLines',idx2,'delimiter','\t');
member.index = cell2mat(memberSection(1));
member.jointA = cell2mat(memberSection(2));
member.jointB = cell2mat(memberSection(3));

frewind(fId);

% Read data from reaction section into structure
reactionSection = textscan(fId, '%d %d %c','headerLines',idx3,'delimiter','\t');
reaction.index = cell2mat(reactionSection(1));
reaction.node = cell2mat(reactionSection(2));
reaction.dir = cell2mat(reactionSection(3));

frewind(fId);

% Read data from external force section into structure
extfSection = textscan(fId, '%d %d %f %c','headerLines',idx4,'delimiter','\t');
extf.index = cell2mat(extfSection(1));
extf.joint = cell2mat(extfSection(2));
extf.force = cell2mat(extfSection(3));
extf.dir = cell2mat(extfSection(4));

frewind(fId);

% Read data from unit section into structure
units = cell2mat(textscan(fId, '%c','headerLines',idx5));
unit = cellstr(units);

fclose(fId) ; % close file

NJ = size(joint.index,1); % number of joints

NM = size(member.index,1); % number of members

NR = size(reaction.index,1); % number of reactions

NE = size(extf.index,1); % number of external forces

if 2*NJ ~= (NM+NR)
    message = '2*(#Joints) does not equal (#Members + #Reactions)';
    fprintf ('\nNJ = %d, NM = %d and NR = %d\n', NJ, NM, NR)
    error(message)
end
% note: add error checking, make sure 2*NJ = NM+NR

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

% PRINT TO FILE
text = '.txt';
prompt = 'Please enter file to write to: ';
fileName = input(prompt,'s');
fileName = strcat(fileName,text);
trust = fopen(fileName,'wt') ;
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
    fprintf(trust,'Member ');
    fprintf(trust,'%02.0f',counter);
    fprintf(trust, ': F= ');
    fprintf(trust, '%10.2f',abs(result(memnum)));
    fprintf(trust, ' ');
    fprintf(trust, units);
    
    if result(memnum)<0
        fprintf(trust, ' [C]');
    else
        fprintf(trust, ' [T]');
    end
    fprintf(trust, '\n');
    counter=counter + 1;
end

for reactnum = 1:NR
    if reaction.dir(reactnum) == 'X'
        fprintf(trust, 'X');
    else
        fprintf(trust, 'Y');
    end
    fprintf(trust,'-direction reaction on joint ');
    fprintf(trust,'%2.0f',reaction.node(reactnum));
    fprintf(trust, ' = ');
    fprintf(trust, '%+10.2f',result(NM+reactnum));
    fprintf(trust, ' ');
    fprintf(trust, units);
    fprintf(trust, '\n');
end

fclose(trust) ; % close file

