#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>

bool inVisited(int which, int a, int *vis)
{
    for(int i=0; i<a; i++)
    {
        if(vis[i]==which)
            return true;
    }
    return false;
}

struct Segment
{
    int from, s;
    int * code;
};

class Entity
{
    int * permutation;

public:
    //numberOfPlaces - number of all places, includes base
    Entity(int numberOfPlaces)
    {
        permutation = new int [numberOfPlaces-1];
        for(int i=0; i<numberOfPlaces-1; i++)
            permutation[i] = 0;

        for(int i=1; i<numberOfPlaces; i++)
        {
            int nr = std::rand()%(numberOfPlaces-1);
            while(permutation[nr]!=0)
            {
                nr++;
                nr = nr%(numberOfPlaces-1);
            }
            permutation[nr] = i;
        }
    }

    void newBorn(int * perm)
    {
        permutation = perm;
    }

    Segment giveCodeSection(int numberOfPlaces)
    {
        Segment s;
        s.from = std::rand()%(numberOfPlaces-2);
        int pom = std::rand()%(numberOfPlaces-2);
        if(pom<s.from)
        {
            s.from = pom;
        }
        s.s = std::rand()%(numberOfPlaces-1-s.from)+1;
        s.code = new int [s.s];
        s.code = &permutation[s.from];
        std::cout<<&permutation[s.from]<<std::endl;
        for(int i=0; i<s.s; i++)
            std::cout<<s.code[i]<<" "<<std::endl;
        return s;
    }

    int * makeUpGenome(int numberOfPlaces, Segment s)
    {
        int * newGenom = new int [numberOfPlaces-1];
        for(int i=0; i<numberOfPlaces-1; i++)
            newGenom[i] = 0;
        for(int i=0; i<s.s; i++)
        {
            std::cout<<s.code[i]<<std::endl;
            newGenom[i+s.from] = s.code[i];
        }
        int j=0;
        int k=0;
        for(int i=0; i<numberOfPlaces-1-s.s; i++)
        {
            while(inVisited(permutation[j],s.s,s.code))
                j++;
            while(s.code[k]==0)
                k++;
            newGenom[k]=permutation[j];
        }
        return newGenom;
    }

    int * retPerm()
    {
        return permutation;
    }

    void swaping(int i, int j)
    {
        int pom = permutation[i];
        permutation[i] = permutation[j];
        permutation[j] = pom;
    }
};

class Evolution
{
    Entity ** generation;
    int numberOfPlaces, population, vehicleCapacity;
    //database is a 2d matrix where the first row is set of distances from BASE
    float ** dataBase;
    Entity * wsk;


public:
    Evolution(int nOP, int vC, int pop = 1000)
    {
        numberOfPlaces = nOP;
        population = pop;
        vehicleCapacity = vC;

        generation = new Entity * [population];
        for(int i=0; i<population; i++)
        {
            wsk = new Entity(numberOfPlaces);
            generation[i] = wsk;
        }

        dataBase = new float * [numberOfPlaces];
        for(int i=0; i<numberOfPlaces; i++)
            dataBase[i] = new float [numberOfPlaces];
    }

    void loadData()
    {
        //poprawiæ, aby poprawnie wczytywa³o tyle danych ile potrzeba
        std::fstream file("distances.csv");
        for(int i=0; i<numberOfPlaces;  i++)
        {
            for(int j=0; j<numberOfPlaces; j++)
            {
                file>>dataBase[i][j];
            }
        }
        file.close();
    }

//////////////////////////
    int * mutation(int * tab)
    {
        int mutated[6];
        for(int i=0 ; i<6; i++)
            mutated[i] = std::rand()%(numberOfPlaces-1);
        for(int i=0; i<3; i++)
        {
            int pom = tab[mutated[i]];
            tab[mutated[i]] = tab[mutated[6-i]];
            tab[mutated[6-i]] = pom;
        }
        return tab;
    }

    void crossOver(int ** couples)
    {
        int j=0;
        for(int i=0; i<population/2; i++)
        {
            while(generation[j]!=NULL)
                j++;
            wsk = new Entity(numberOfPlaces-1);
            generation[j] = wsk;
            wsk->newBorn(mutation(generation[couples[i/2][i%2]]->makeUpGenome(numberOfPlaces,generation[couples[i/2][(i+1)/2]]->giveCodeSection(numberOfPlaces))));
        }
    }
////////////////////////
    int ** mixing(bool * listOfKilled)
    {
        int numberOfCouples = population/4;
        int ** couples = new int * [numberOfCouples];
        for(int i=0; i<numberOfCouples; i++)
        {
            couples[i] = new int [2];
            int x=i;
            while(listOfKilled[x])
                x++;
            listOfKilled[x] = true;
            couples[i][0]=x;
            x++;
            int luck = std::rand()%5;
            while((listOfKilled[x]==true)||(luck!=0))
            {
                luck = std::rand()%5;
                x++;
                if(x==population)
                    x=couples[i][0]+1;
            }
            listOfKilled[x] = true;
            couples[i][1]=x;
        }

        return couples;
    }

    bool * killing()
    {
        int numberOfKilled = population/2;
        bool * listOfKilled = new bool [population];
        for(int i=0; i<population; i++)
            listOfKilled[i] = false;

        for(int i=0; i<numberOfKilled; i++)
        {
            int j=0;
            int luck = std::rand()%(population/7);
            while((listOfKilled[population-1-j]==true)||(luck!=0))
            {
                luck = std::rand()%(population/7);
                j++;
                j=j%(population);
            }
            listOfKilled[population-1-j] = true;
            delete generation[population-1-j];
            generation[population-1-j] = NULL;
        }
        return listOfKilled;
    }
/////////////////////////////////////////////
    void swaping(int i, int j, Entity * e)
    {
        float distance = 0;
        if((i+1)%vehicleCapacity==0)
        {
            distance+=dataBase[e->retPerm()[i-1]][e->retPerm()[j]];
            distance+=dataBase[e->retPerm()[j]][0];
            distance-=dataBase[e->retPerm()[i-1]][e->retPerm()[i]];
            distance-=dataBase[e->retPerm()[i]][0];
        }
        else if((i+1)%vehicleCapacity==1)
        {
            distance+=dataBase[0][e->retPerm()[j]];
            distance+=dataBase[e->retPerm()[j]][e->retPerm()[i+1]];
            distance-=dataBase[0][e->retPerm()[i]];
            distance-=dataBase[e->retPerm()[i]][e->retPerm()[i+1]];
        }
        else
        {
            distance+=dataBase[e->retPerm()[i-1]][e->retPerm()[j]];
            distance+=dataBase[e->retPerm()[j]][e->retPerm()[i+1]];
            distance-=dataBase[e->retPerm()[i-1]][e->retPerm()[i]];
            distance-=dataBase[e->retPerm()[i]][e->retPerm()[i+1]];
        }

        if(j+2==numberOfPlaces)
        {
            distance+=dataBase[e->retPerm()[j-1]][e->retPerm()[i]];
            distance+=dataBase[e->retPerm()[i]][e->retPerm()[0]];
            distance-=dataBase[e->retPerm()[j-1]][e->retPerm()[j]];
            distance-=dataBase[e->retPerm()[j]][e->retPerm()[0]];
        }
        else
        {
            if((j+1)%vehicleCapacity==0)
            {
                distance+=dataBase[e->retPerm()[j-1]][e->retPerm()[i]];
                distance+=dataBase[e->retPerm()[i]][0];
                distance-=dataBase[e->retPerm()[j-1]][e->retPerm()[j]];
                distance-=dataBase[e->retPerm()[j]][0];
            }
            else if((j+1)%vehicleCapacity==1)
            {
                distance+=dataBase[0][e->retPerm()[i]];
                distance+=dataBase[e->retPerm()[i]][e->retPerm()[j+1]];
                distance-=dataBase[0][e->retPerm()[j]];
                distance-=dataBase[e->retPerm()[j]][e->retPerm()[j+1]];
            }
            else
            {
                distance+=dataBase[e->retPerm()[j-1]][e->retPerm()[i]];
                distance+=dataBase[e->retPerm()[i]][e->retPerm()[j+1]];
                distance-=dataBase[e->retPerm()[j-1]][e->retPerm()[j]];
                distance-=dataBase[e->retPerm()[j]][e->retPerm()[j+1]];
            }
        }

        if(distance<0)
        {
            e->swaping(i,j);
        }

    }

    float computeDistance(Entity e)
    {
        float distance = 0;
        distance += dataBase[0][e.retPerm()[0]];
        distance += dataBase[numberOfPlaces-1][0];
        for(int i=1; i<numberOfPlaces-1; i++)
        {
            if(i%vehicleCapacity==0)
            {
                distance+=dataBase[e.retPerm()[i-1]][0];
                distance+=dataBase[0][e.retPerm()[i]];
            }
            else
                distance+=dataBase[e.retPerm()[i-1]][e.retPerm()[i]];
        }
        return distance;
    }

    void twoOpt(Entity * e)
    {
        for(int i=0; i<numberOfPlaces-2; i++)
        {
            for(int j=i+1; j<numberOfPlaces-1; j++)
            {
                swaping(i,j,e);
            }
        }
    }

    void scal(int k, int kmax, int l, int lmax, Entity ** gen, float * score, float * newScore)
    {
        int i=k;
        while((k<=kmax)&&(l<=lmax))
        {
            if(score[k]<=score[l])
            {
                gen[i] = generation[k];
                newScore[i] = score[k];
                k++;
            }
            else
            {
                gen[i] = generation[l];
                newScore[i] = score[l];
                l++;
            }
            i++;
        }

        int m;
        if(k>kmax)
            m=l;
        else
            m=k;

        while(i<=lmax)
        {
            gen[i] = generation[m];
            newScore[i] = score[m];
            m++;
            i++;
        }
    }

    void mergeSort(float * score)
    {
        for(int i=0; i<ceil(log(population))+2; i++)
        {
            float * newScore = new float [population];
            Entity ** sortedGeneration = new Entity * [population];

            for(int j=0; j<ceil(population/pow(2,i+1)); j++)
            {
                int k = j*pow(2,i+1);
                int l = (j+1)*pow(2,i+1);
                if(l>=population)
                    continue;
                int kmax = l-1;
                int lmax = l+pow(2,i+1)-1;
                if(lmax>=population)
                    lmax = population-1;

                scal(k,kmax,l,lmax,sortedGeneration,score,newScore);
            }
            score = newScore;
            generation = sortedGeneration;
        }
        //for(int j=0; j<population; j++)
            //std::cout<<j<<") "<<computeDistance(*generation[j])<<std::endl;
    }

    void calcAndSort()
    {
        float * score = new float [population];
        for(int i=0; i<population; i++)
        {
            twoOpt(generation[i]);

            score[i] = computeDistance(*generation[i]);
        }
        mergeSort(score);
    }

    void pushEvolution(int times=1)
    {
        for(int i=0; i<times; i++)
        {
            calcAndSort();
            crossOver(mixing(killing()));
        }
        calcAndSort();
        std::cout<<" the best: "<<computeDistance(*generation[0]);
    }
};

int main()
{
    srand( time( NULL ) );
    Evolution evo(44,8);
    evo.loadData();
    evo.pushEvolution(1);
    return 0;
}
