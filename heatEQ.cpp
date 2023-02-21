/* This file contains several implementations of the heat equation using explicit and implicit methods, some of them 
work with the standard vectors and others with arrays. The goal was to have a refresher on the basics of scientific computing and 
help me start to learn how to code in c++ */

#include<vector>
#include<cmath>
#include<iostream>

// set initial distribution
float init_dist (float x){
    return std::pow(x,2);
};


// set the forcing function
float forcing(float x, float t){
    return 0.0;
};


// This section defines several useful linear algebra operations that work on the standard library vectors

/* function to swap either two given columns of rows within a matrix. Swap is done in place. If row=True it will swap rows a and b. 
If row=False then the columns will be swapped */
void matrix_row_swap(std::vector<std::vector<float>> &A, bool row, int a, int b, int l ){
    float dum;
    if (row){
        for (int i=0; i<l; i++){
            dum = A[a][i];
            A[a][i] = A[b][i];
            A[b][i] = dum;
        };
    }else{
        for(int i=0; i<l; i++){
            dum = A[i][a];
            A[i][a] = A[i][b];
            A[i][b] = dum;
        };
    };
};


// take a vector and transform it in place by a given matrix
void vector_map (std::vector<float> &v, std::vector<std::vector<float>> &A, int l ){
    std::vector<float> next;
    float val;

    // compute the next step of the solution
    next.push_back(v[0]);
    for (int i=1; i<l-1; i++){
        val = 0;
        for(int j=0; j<l; j++){
            val += A[i][j] * v[i];
            };
        next.push_back(val);
    };
    next.push_back(v[l-1]);

    // update the distribution
    for (int i=0; i<l; i++){v[i] = next[i];};

};


// function to take a matrix and compute the inverse in place
void gaussj (std::vector<std::vector<float>> &A, int l ){

    float inv;
    float reduce;

    // augment the matrix with the identity
    for(int i=0; i<l; i++){
        for(int j=0; j<l; j++){
            if (i == j) {
                 A[i].push_back(1);
            }else{ 
                A[i].push_back(0);
            };
        };
    };

    // loop through and create the inverse matrix
    for (int i=0; i<l; i++){

        // normalize the current row of the matrix
        inv = 1.0/A[i][i];
        for (int j=0; j<2*l; j++){
            A[i][j] *= inv;
        };

        // reduce the other rows of the matrix
        for (int j=0; j<l; j++){
            if (j!=i){
                reduce = A[j][i];
                for(int k=0; k<2*l; k++){
                    A[j][k] -= A[i][k]*reduce;
                };
            };
        };
    };

    // delete the augmented part of the matrix and copy over the inverse
    for (int i=0; i<l; i++){
        for (int j=l-1; j>=0; j--){
            A[i][j] = A[i][j+l];
            A[i].pop_back();
        };
    };

}; 


/* Below are several classes that may be used to simulate heat transfer in one dimension. Please note that only the 
fast_heatEQ class accepts a source/forcing function and all simulations use dirichlet boundary conditions */


/* class that is used to simulate heat transfer. The child classes correspond to the various numerical methods that may 
be used to run the simulation */
class heateq{

    public:

        // prints out the current solution values and time
        void display(){
            for (float i : heat_dist){
                std::cout << i << ' ' ;
            };
            std::cout << '\n';
        };

        // then initialize the distribution and parameters acordingly
        heateq(int num_points_, float time_, float dt_, float (*dist)(float), 
        float left_endpoint_, float right_endpoint_, float left_condition, float right_condition ){

            reset_sim(num_points_, time_, dt_, dist, left_endpoint_, 
            right_endpoint_, left_condition, right_condition );
        };

        // reset the time, boundary conditions and heat distribution
        void reset_sim (int num_points_, float time_, float dt_, float (*dist)(float), 
        float left_endpoint_, float right_endpoint_, float left_condition_, float right_condition_ ) {

           // set parameters
            num_points = num_points_;
            time = time_;
            dt = dt_;
            left_endpoint = left_endpoint_;
            right_endpoint = right_endpoint_;
            left_condition = left_condition_;
            right_condition = right_condition_;
            dx = (right_endpoint - left_endpoint)/(num_points-1.0);

            // initial distribution
            heat_dist.push_back(left_condition);
            for(int i=1; i<num_points-1; i++){
                heat_dist.push_back(dist(left_endpoint + dx*(i*1.0)));
            };
            heat_dist.push_back(right_condition);
            display();
        };

    protected:
        int num_points;
        float left_endpoint;
        float right_endpoint;
        float left_condition;
        float right_condition;
        float time;
        float dt;
        float dx;

        // the current distribution on the line
        std::vector<float> heat_dist; 

};


// childs class that uses an explicit finite difference method
class explct : public heateq{

    public:

        // class constructor
        explct (int num_points_, float time_, float dt_, float (*dist)(float), 
            float left_endpoint_, float right_endpoint_, float left_condition, float right_condition)
            : heateq(num_points_, time_, dt_, dist, left_endpoint_, 
            right_endpoint_, left_condition, right_condition ) {};

        // run a step of the simulation
        void next_step (){
        
            std::vector<float> next;

            // handle endpoint case
            next.push_back(left_condition);
            // main loop for intermediate values
            for (int i=1; i<num_points -1; i++ ){
                next.push_back(heat_dist[i] + dt*(heat_dist[i-1] - 2.0*heat_dist[i] + heat_dist[i+1])/std::pow(dx,2));
            };
            // handle endpoint case
            next.push_back(right_condition);

            // copy values to the distribution
            for (int i=0; i<num_points; i++){
                heat_dist[i] = next[i];
            };
            time++;
        };

};


// child class that uses an implicit finite difference method
class implct : public heateq{

    public:

        implct (int num_points_, float time_, float dt_, float (*dist)(float), 
            float left_endpoint_, float right_endpoint_, float left_condition, float right_condition)
            : heateq(num_points_, time_, dt_, dist, left_endpoint_, 
            right_endpoint_, left_condition, right_condition ) {

                // use gauss joran matrix inversion to determine the step matrix
                float lmbda = dt/std::pow(dx,2);
                for(int i=0; i<num_points; i++){
                    std::vector<float> next_row;
                    
                    for(int j=0; j<num_points; j++){
                        if (j == i-1){next_row.push_back(-lmbda);}
                        else if(j == i){next_row.push_back(1 + lmbda*2.0);}
                        else if(j == i+1){next_row.push_back(-lmbda);}
                        else{next_row.push_back(0);};  
                    };
                    step_matrix.push_back(next_row);
                };
                gaussj(step_matrix, num_points);
            };

        void next_step (){
            vector_map(heat_dist, step_matrix, num_points);
        };

        // print out the step matrix
        void disp_matrix(){
            for(int i=0; i<num_points; i++){
                for(int j=0; j<num_points; j++){
                    std::cout<<step_matrix[i][j]<<' ';
                };
                std::cout<<'\n';
            };
            std::cout<<'\n';
        };


    private:
        std::vector<std::vector<float>> step_matrix;
};


// child class that uses the Crank-Nicholson
class crank_nicholson : public heateq{

    public:

        crank_nicholson (int num_points_, float time_, float dt_, float (*dist)(float), 
            float left_endpoint_, float right_endpoint_, float left_condition, float right_condition)
            : heateq(num_points_, time_, dt_, dist, left_endpoint_, 
            right_endpoint_, left_condition, right_condition ) {


                // use gauss joran matrix inversion to determine the step matrix
                float lmbda = dt/std::pow(dx,2);
                for(int i=0; i<num_points; i++){
                    std::vector<float> next_row;
                    
                    for(int j=0; j<num_points; j++){
                        if (j == i-1){next_row.push_back(-lmbda);}
                        else if(j == i){next_row.push_back(2 + lmbda*2.0);}
                        else if(j == i+1){next_row.push_back(-lmbda);}
                        else{next_row.push_back(0);};  
                    };
                    step_matrix.push_back(next_row);
                };
                gaussj(step_matrix, num_points);
                
                //initialize the stepping vector
                step_vector.push_back(left_condition);
                for (int j=1; j<num_points-1; j++){ step_vector.push_back(0.0);};
                step_vector.push_back(right_condition);

            };


        void next_step (){
            // update the step vector
            float lmbda = dt/std::pow(dx,2);
            for (int i=1; i<num_points-1; i++){
                step_vector[i] = lmbda*heat_dist[i-1] + (2 - 2*lmbda)*heat_dist[i] + lmbda*heat_dist[i+1];
            };
            //multiply to get next step of solution
            vector_map(step_vector, step_matrix, num_points);
            for (int i=1; i<num_points-1; i++){heat_dist[i] = step_vector[i];};
        };

        // print out the step matrix
        void disp_matrix(){
            for(int i=0; i<num_points; i++){
                for(int j=0; j<num_points; j++){
                    std::cout<<step_matrix[i][j]<<' ';
                };
                std::cout<<'\n';
            };
            std::cout<<'\n';
        };

    private:
        std::vector<std::vector<float>> step_matrix;
        std::vector<float> step_vector;
};


// a much more optimized implimentation of the heat equation, uses crank nicholson method for high accuracy without too much work.
class fast_heateq{

    public:
        fast_heateq(float (*dist)(float), float (*forcing)(float, float), float left_condition_, float right_condition_,
         float left_endpoint_, float right_endpoint_, float dt_, float dx_){

            // initialize parameters
            left_condition = left_condition_;
            right_condition = right_condition_;
            left_endpoint = left_endpoint_;
            right_endpoint = right_endpoint_;
            dt = dt_;
            dx = dx_;
            t = 0.0;

            num_points = (right_endpoint - left_endpoint) / dx + 1;
            lmbda = dt/std::pow(dx,2);
            a = 2.0 + 2.0*lmbda;
            b = -lmbda;
            gamma = b/a;

            // initialize distribution
            heat_dist = new float[num_points];
            for (int i=1; i<num_points-1; i++){
                heat_dist[i] = dist(left_endpoint + dx*i);
            };
            heat_dist[0] = left_condition;
            heat_dist[num_points-1] = right_condition;
        };


        // compute the solution at the next time step
        void next_step(){

            // update heat distribution to the weighted sums used for the scheme
            float uprev = heat_dist[0];
            float ucur = heat_dist[1];
            float unext = heat_dist[2];
            float beta;

            heat_dist[0] = left_condition;
            for (int i=1; i< num_points-1; i++){
                heat_dist[i] = lmbda*uprev + (2.0 - 2.0*lmbda)*ucur + lmbda*unext + 0.5*(forcing(heat_dist[i], t) + forcing(heat_dist[i], t+dt));
                uprev = ucur;
                ucur = unext;
                unext = heat_dist[i+2];
            };
            heat_dist[num_points-1] = right_condition;

            // compute the foreward substitution
            for (int i=1; i<num_points-1; i++){
                beta = a - b*gamma;
                heat_dist[i] = (heat_dist[i] - b*heat_dist[i-1])/a;
            };

            // compute the backward substitution
            for (int i=num_points-2; i>0; i--){
                heat_dist[i] -= gamma*heat_dist[i+1];
            };

            // update simulation time
            t+=dt;
        };


        //print out the current heat distribution
        void display(){
            std::cout << "(";
            for (int i=0; i<num_points; i++){
                std::cout<<heat_dist[i]<<", ";
            };
            std::cout << "), \n";
        };


        //class destructor, eliminates the dynamically allocated arrays
        ~ fast_heateq(){
            delete [] heat_dist;
        };


    private:
    // simulation parameters
    float dx;
    int num_points;
    float left_endpoint;
    float right_endpoint;
    float left_condition;
    float right_condition;
    float dt;
    float lmbda;
    float a;
    float b;
    float gamma;
    float t;

    // pointer to the array used for computation
    float* heat_dist;
    
};



// run some basic checks on the classes
int main(){

    // precision of discretization
    float dt = 0.00005;
    float dx = 0.01;

    // window 
    float left_endpoint = 0.0;
    float right_endpoint = 1.0;

    // boundary conditions (dirichlet)
    float left_condition = 0.0;
    float right_condition = 0.0;


    fast_heateq test(init_dist, forcing, left_condition, right_condition,
                     left_endpoint, right_endpoint, dt, dx);

    for (int i=0; i<100; i++){
        if (i%20 == 0){
            test.display();
        };
        test.next_step();
    };

}
