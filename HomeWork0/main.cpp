 #include<cmath>
 #include<iostream>
 #include<eigen3/Eigen/Core>
 #include<eigen3/Eigen/Dense>

int main(){
    std::cout << "HomeWork \n";

    Eigen::Vector3f p(2.0f,1.0f,1.0f);
    Eigen::Matrix3f r45Move;
    r45Move << std::cos(45.0/180.0*acos(-1)),std::sin(-45.0/180.0*acos(-1)),1.0f,
    std::sin(45.0/180.0*acos(-1)),std::cos(45.0/180.0*acos(-1)),2.0f,
    0.0f,0.0f,1.0f;
        
    std::cout << r45Move * p << std::endl;
}