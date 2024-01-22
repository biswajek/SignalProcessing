pkg load statistics;
pkg load communications;
close('all');
clear;

% Generate and plot dataset X
randn('seed',0);
m1=[1, 1]'; m2=[3, 3]';m3=[2, 6]';
m=[m1 m2 m3];

S(:,:,1)=0.1*eye(2);
S(:,:,2)=0.2*eye(2);
S(:,:,3)=0.3*eye(2);
P=[0.4 0.4 0.2];
N=1000;
sed=0;
rand('seed',sed);
[l,c]=size(m);

%Constructing the P_acc vector. This is necessary for picking randomly
%one of the c normal distributions in order to generate a point.
P_acc=P(1);
for i=2:c
    t=P_acc(i-1)+P(i);
    P_acc=[P_acc t];
end

% Generation of the data set
X=[];
y=[];
for i=1:N
    t=rand;
    ind=sum(t>P_acc)+1;  % Index of the normal distribution that will generate the i-th vector
    X=[X; mvnrnd(m(:,ind)',S(:,:,ind),1)];
    y=[y ind];
end
X=X';


[l,N]=size(X); % N=no. of data vectors, l=dimensionality
[l,c]=size(m); % c=no. of classes
figure(1);
hold on;
for i=1:N
    plot(X(1,i),X(2,i),'.r')
    hold on
end

% Plot of the class centroids
##for j=1:c
##    plot(m(1,j),m(2,j),'*r',"markersize",20,"LineWidth",1.5)
##    hold on
##end


m1_ini=[0; 2];m2_ini=[5; 2];m3_ini=[5; 5];
m=[m1_ini m2_ini m3_ini];

% Plot of the class centroids
for j=1:c
    scatter(m(1,j),m(2,j),'m',"filled")
end

s=[.15 .27 .4];
Pa=[1/3 1/3 1/3];
e_min=10^(-5);

x=X';
m=m';
[p,n]=size(x);
[J,n]=size(m);

e=e_min+1;

Q_tot=[];
e_tot=[];

iter=0;
while (e>e_min)
    iter=iter+1;
    e;

    P_old=Pa;
    m_old=m;
    s_old=s;

    % Determine P(j|x_k; theta(t))
    for k=1:p
        temp=gauss(x(k,:),m,s);
        P_tot=temp*Pa';
        for j=1:J
            P(j,k)=temp(j)*Pa(j)/P_tot;
        end
    end

    % Determine the log-likelihood
    Q=0;
    for k=1:p
        for j=1:J
            Q=Q+P(j,k)*(-(n/2)*log(2*pi*s(j)) - sum( (x(k,:)-m(j,:)).^2)/(2*s(j)) + log(Pa(j)) );
        end
    end
    Q_tot=[Q_tot Q];

    % Determine the means
    for j=1:J
        a=zeros(1,n);
        for k=1:p
            a=a+P(j,k)*x(k,:);
        end
        m(j,:)=a/sum(P(j,:));
    end
    mt=m';
    figure(1);
    [l,c]=size(mt); % c=no. of classes
    for j=1:c
      plot(mt(1,j)',mt(2,j)','+m',"LineWidth",1.5)
      %hold on
    end
    hold on;
    %scatter(m(1,:)',m(2,:)','k',"filled");
    % Determine the variances
    for j=1:J
        b=0;
        for k=1:p
            b=b+ P(j,k)*((x(k,:)-m(j,:))*(x(k,:)-m(j,:))');
        end
        s(j)=b/(n*sum(P(j,:)));

        if(s(j)<10^(-10))
            s(j)=0.001;
        end

    end

    % Determine the a priori probabilities
    for j=1:J
        a=0;
        for k=1:p
            a=a+P(j,k);
        end
        Pa(j)=a/p;
    end

    e=sum(abs(Pa-P_old))+sum(sum(abs(m-m_old)))+sum(abs(s-s_old));
    e_tot=[e_tot e];
end
m=m';
figure(1);
[l,c]=size(m); % c=no. of classes
for j=1:c
    plot(m(1,j)',m(2,j)','+b',"markersize",12,"LineWidth",2)
    hold on
end
