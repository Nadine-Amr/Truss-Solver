clear
clc
n=input('Please enter the number of nodes:');
m=input('Please enter the number of members:');
if m>(3*n)-6
    fprintf('The structure is indeterminte.');
else
    initialfolder= 'C:\Users\acer';
    folderchoice=uigetfile(initialfolder);
    narray=xlsread(folderchoice,1);
    marray=xlsread(folderchoice,2);
    extforces=xlsread(folderchoice,3);
    rarray=xlsread(folderchoice,4);
    if sum(sum(rarray(:,2:4)))>6
        fprintf('The structure is indeterminte.');
    else
        [k,c]=size(extforces);
        fex=sym(zeros(n,3));
        [s,d]=size(rarray);
        counter=0;
        for i=1:s
            for j=2:d
                if rarray(i,j)==1
                    counter=counter+1;
                    syms (strcat('r',num2str(counter)))
                    fex(rarray(i,1),(j-1))=strcat('r',num2str(counter));
                end
            end
        end
        for i=1:k
            fex(extforces(i,1),1)=extforces(i,2)*cosd(extforces(i,3))*cosd(extforces(i,4));
            fex(extforces(i,1),2)=extforces(i,2)*cosd(extforces(i,3))*sind(extforces(i,4));
            fex(extforces(i,1),3)=extforces(i,2)*sind(extforces(i,3));
        end
        summation=sum(fex);
        eqn1=summation(1)==0;
        eqn2=summation(2)==0;
        eqn3=summation(3)==0;
        eqn4=0;
        eqn5=0;
        eqn6=0;
        for i=1:n
            eqn4=eqn4+((fex(i,3))*(narray(i,2)-narray(1,2)))-((fex(i,2))*(narray(i,3)-narray(1,3)));
            eqn5=eqn5-((fex(i,3))*(narray(i,1)-narray(1,1)))+((fex(i,1))*(narray(i,3)-narray(1,3)));
            eqn6=eqn6+((fex(i,2))*(narray(i,1)-narray(1,1)))-((fex(i,1))*(narray(i,2)-narray(1,2)));
        end
        eqn4=eqn4==0;
        eqn5=eqn5==0;
        eqn6=eqn6==0;
        reactions=solve(eqn1,eqn2,eqn3,eqn4,eqn5,eqn6);
        reactions=struct2cell(reactions);
        reactionsarr=zeros(length(reactions),1);
        for i=1:length(reactions)
            reactionsarr(i)=reactions{i};
        end
        counter=0;
        for i=1:s
            for j=2:d
                if rarray(i,j)==1
                    counter=counter+1;
                    fex(rarray(i,1),(j-1))=reactionsarr(counter);
                end
            end
        end
        for i=1:m
            syms (strcat('F',num2str(i)))
        end
        syms eqn
        numeqns=0;
        eqn=sym(zeros(1,1));
        for i=1:n       
            marray2=(marray(:,2:3))';
            z1= find(sum(marray2==i));
            tonodes=sum(marray2(1:2,z1))-i;
            A=zeros(3,length(tonodes));
            syms x
            x=sym(zeros(1,length(tonodes)));
            for j=1:length(tonodes)
                len=sqrt((narray(tonodes(j),1)-narray(i,1))^2+(narray(tonodes(j),2)-narray(i,2))^2+(narray(tonodes(j),3)-narray(i,3))^2);
                A(1,j)=(narray(tonodes(j),1)-narray(i,1))/len;
                A(2,j)=(narray(tonodes(j),2)-narray(i,2))/len;
                A(3,j)=(narray(tonodes(j),3)-narray(i,3))/len;
                x(j)=strcat('F',num2str(z1(j)));
            end
            B=-fex(i,:)';
            for j=1:3
                numeqns=numeqns+1;
                eqn(numeqns,1)=0;
                for k=1:length(tonodes)
                    eqn(numeqns,1)=eqn(numeqns,1)+((A(j,k))*x(k));
                end
                eqn(numeqns,1)=eqn(numeqns,1)==B(j);
            end
        end
        fin=solve(eqn(:,1));
        fin=struct2cell(fin);
        intforces=zeros(length(fin),1);
        for i=1:length(fin)
            intforces(i)=fin{i};
        end
        intforces
    end
end
