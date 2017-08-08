#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include <mpi.h>
#include <string.h>

struct Term{
    double Cofficient;
    double Power;
};

int numOfTerms;
struct Term Equation[100];
void swap(double *x,double *y){
    double tmp=*x;
    *x=*y;
    *y=tmp;
}
int isDigit(char c)
{
	return c >= '0'&&c <= '9';
}

int isValid(char equation[]){
	int sz = strlen(equation), i,point=0,power=0;
	equation[sz-1]='\0';
	sz--;
	if(sz==0)
        return 0;
	for (i = 0; i < sz; ++i){
	if (equation[i] == ' '){
			memmove(&equation[i], &equation[i + 1], strlen(equation) - i);
			sz--;
			i--;
		}
		else if (!isDigit(equation[i]) && equation[i] != 'X'&&equation[i] != '+'&&equation[i] != '-'&&equation[i] != '^'&&equation[i] != '.')
			return 0;
		if (i == sz - 1 && !isDigit(equation[i])&& equation[i] != 'X')
			return 0;
		if (equation[i] == '-'&&equation[i + 1] == '+'){
			memmove(&equation[i+1], &equation[i+2], strlen(equation) - i-1);
			sz--;
			i--;
		}
		else if (equation[i] == '+'&&equation[i + 1] == '-'){
			memmove(&equation[i], &equation[i + 1], strlen(equation) - i );
			sz--;
			i--;
		}
		else if (equation[i] == '+'&&equation[i + 1] == '+'){
			memmove(&equation[i], &equation[i + 1], strlen(equation) - i);
			sz--;
			i--;
		}
		else if (equation[i] == '-'&&equation[i + 1] == '-'){
			memmove(&equation[i + 1], &equation[i + 2], strlen(equation) - i - 1);
			equation[i] = '+';
			sz--;
			i--;
		}

		if (equation[i] == '.'){
            if(i==0||(i>0&&!isDigit(equation[i-1]))||i==sz-1||(i<sz-1&&!isDigit(equation[i+1])))
                return 0;
			if (!point)
				point = 1;
			else
				return 0;
		}
		else if (!isDigit(equation[i]))
			point = 0;

		if (equation[i] == '^'){
            if(i==0||(i>0&&equation[i-1]!='X')||i==sz-1||(i<sz-1&&!isDigit(equation[i+1])))
                return 0;
			if (!power)
				power = 1;
			else
				return 0;
		}
		else if (!isDigit(equation[i]))
			power = 0;
	}
	return 1;
}

void extractTerms(char equation[]){
	char term[100];
	int sz = strlen(equation), i,k;
	for (i = 0,numOfTerms=0,k=0; i<sz; ++i){
		if (isDigit(equation[i])||equation[i]=='.'){
			term[k++] = equation[i];
			term[k] = '\0';
		}
		else if (equation[i] == 'X'){
			if (k == 0)
			{
				term[k++] = '1';
				term[k] = '\0';
			}
			if (Equation[numOfTerms].Cofficient==-1)
				Equation[numOfTerms].Cofficient = -atof(term);
			else
				Equation[numOfTerms].Cofficient= atof(term);

			Equation[numOfTerms].Power = 1;
			memset(term, 0, sizeof(term));
			k = 0;
			term[k] = '\0';
		}
		else if (equation[i] == '+'){
			if (k == 0)
			{
				if (i != 0)
					numOfTerms++;
				continue;
			}
			if (Equation[numOfTerms].Power == -1)
				Equation[numOfTerms].Power = -atof(term);
			else
				Equation[numOfTerms].Power = atof(term);
			memset(term, 0, sizeof(term));
			k = 0;
			if (i != 0)
				numOfTerms++;
			term[k] = '\0';
		}
		else if (equation[i] == '-'){
			if (k == 0)
			{
				if (i != 0)
					numOfTerms++;
				Equation[numOfTerms].Cofficient = -1;
				continue;
			}
			if (Equation[numOfTerms].Power == -1)
				Equation[numOfTerms].Power = -atof(term);
			else
				Equation[numOfTerms].Power = atof(term);
			memset(term, 0, sizeof(term));
			k = 0;
			if (i != 0)
				numOfTerms++;
			Equation[numOfTerms].Cofficient = -1;
			term[k] = '\0';
		}
	}
	if (strlen(term)>0)
	if (Equation[numOfTerms].Cofficient == -1)
		Equation[numOfTerms++].Cofficient = -atof(term);
	else if(Equation[numOfTerms].Cofficient == 0)
		Equation[numOfTerms++].Cofficient = atof(term);
    else
		Equation[numOfTerms++].Power = atof(term);


}

double f(double x){
	double result = 0;
	int i;
	for (i = 0; i<=numOfTerms; ++i){
		result += Equation[i].Cofficient * pow(x, Equation[i].Power);
	}
	return result;
}

double min(double x, double y){
	return x>y ? y : x;
}



int main(int argc, char** argv){

    char equation[100];
    MPI_Init(&argc,&argv);
    struct Term _Term;
    int count;
    count = 1;
    MPI_Datatype array_of_types[count];
    array_of_types[0] = MPI_DOUBLE;

    int array_of_blocklengths[count];
    array_of_blocklengths[0] = 2;

    MPI_Aint array_of_displaysments[count];
    MPI_Aint address1, address2;
    MPI_Get_address(&_Term,&address1);
    MPI_Get_address(&_Term.Cofficient,&address2);
    array_of_displaysments[0] = address2 - address1;

    /*Create MPI Datatype and commit*/
    MPI_Datatype MPI_TERM;
    MPI_Type_create_struct(count, array_of_blocklengths, array_of_displaysments, array_of_types, &MPI_TERM);
    MPI_Type_commit(&MPI_TERM);

    double
        low,
        high,
        dx,
        x,
        y,
        yold,
        ans=0,msg=0,diff;

    int rank,i,size,n;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if(rank==0){
        system("clear");
        printf("Welcome to Integration Calculator\n\nF(X) = ");
        fgets(equation, 100, stdin);
        printf("Please enter the lower limit : ");
        scanf("%lf",&low);
        printf("Please enter the upper limit : ");
        scanf("%lf",&high);
        printf("Please enter the number of trapezoides ( Control the accuracy ) : ");
        scanf("%d",&n);
        if(low>high)
            swap(&low,&high);
        if(isValid(equation)==0)
        {
            puts("Invalid Input !");
            return 0;
        }
        extractTerms(equation);
        dx=(high-low)/n;
        if(dx==0)
            dx++;
    }
    MPI_Bcast(&low,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&high,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&dx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Equation,100,MPI_TERM,0,MPI_COMM_WORLD);
    MPI_Bcast(&numOfTerms,1,MPI_INT,0,MPI_COMM_WORLD);
    if(rank==0){
        diff=high-low;
        low=low+rank*diff/size;
        high=low+diff/size;
        yold=f(low);
        for(x=low+dx;x<=high;x+=dx){
            y=f(x);
            ans+=min(yold,y)*dx+((dx)/2.0*fabs(yold-y));
            yold=y;
        }
        x-=dx;
        dx=high-x;
        y=f(x+dx);
        ans+=min(yold,y)*dx+((dx)/2.0*fabs(yold-y));
        for(i=1;i<size;++i){
            MPI_Recv(&msg,1,MPI_DOUBLE,MPI_ANY_SOURCE,5,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
            ans+=msg;
        }
       printf("Answer : %f\n",ans);
    }
    else {
        diff=high-low;
        low=low+rank*diff/size;
        high=low+diff/size;
        yold = f(low);
        for(x=low+dx;x<=high;x+=dx){
            y=f(x);
            ans+=min(yold,y)*dx+((dx)/2.0*fabs(yold-y));
            yold=y;
        }
        x-=dx;
        dx=high-x;
        y=f(x+dx);
        ans+=min(yold,y)*dx+((dx)/2.0*fabs(yold-y));
        MPI_Send(&ans,1,MPI_DOUBLE,0,5,MPI_COMM_WORLD);

    }

    MPI_Finalize();
	return 0;
}
