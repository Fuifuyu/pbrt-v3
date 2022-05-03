#pragma once

#ifndef MY_TYPES_H
#define MY_TYPES_H

#include <array>
#include <vector>
#include <unordered_map>
#include <mutex>

namespace MyTypes{
    template<size_t N>
    class arrayd{
        std::array<double, N> array{0};
    public:
        arrayd(){}
        arrayd(std::array<double, N>  a) : array(a){}
        arrayd(double x, double y) {
            array[0] = x;
            array[1] = y;
        }
        
        double &x() {return array[0];}
        double const &x() const {return array[0];}
        double &y() {return array[1];}
        double const &y() const {return array[1];}
        double &z() {
            if(N<3) 0;
            return array[2];
        }
        double const &z() const {
            if(N<3) 0;
            return array[2];
        }
        int n() const {
            return N;
        }
        
        double &operator[](int i) {return array[i];}
        double const &operator[](int i) const {return array[i];};
        
        arrayd operator-(arrayd const &b) const {
            if(N != b.n()) throw exception("Operation between different dimension.");
            arrayd<N> res(array);
            if(N == 2) {
                res[0] -= b[0]; res[1] -= b[1];
                return res;
            }
            if(N == 3) {
                res[0] -= b[0]; res[1] -= b[1]; res[2] -= b[2];
                return res;
            }
            for(int i = 0; i < N; i++){
                res[i] -= b[i];
            }
            return res;
        }
        
        arrayd operator+(arrayd const &b) const{
            if(N != b.n()) throw exception("Operation between different dimension.");
            arrayd<N> res(array);
            if(N == 2) {
                res[0] += b[0]; res[1] += b[1];
                return res;
            }
            if(N == 3) {
                res[0] += b[0]; res[1] += b[1]; res[2] += b[2];
                return res;
            }
            for(int i = 0; i < N; i++){
                res[i] += b[i];
            }
            return res;
        }
        
        arrayd& operator+=(arrayd<N> const &a){
            if(N==2){
                this->x() += a[0]; this->y() += a[1];
                return *this;
            }
            if(N==3){
                this->x() += a[0]; this->y() += a[1]; this->z() += a[2];
                return *this;
            }
            throw std::exception("Not implemented");
        }

        arrayd& operator*=(double a){
            if(N==2){
                this->x() *= a; this->y() *= a;
                return *this;
            }
            if(N==3){
                this->x() *= a; this->y() *= a; this->z() *= a;
                return *this;
            }
            throw std::exception("Not implemented");
        }
        arrayd operator*(double b) const {
            arrayd<N> res(array);
            if(N == 2) {
                res[0] *= b; res[1] *= b;
                return res;
            }
            if(N == 3) {
                res[0] *= b; res[1] *= b; res[2] *= b;
                return res;
            }
            for(int i = 0; i < N; i++){
                res[i] *= b;
            }
            return res;
        }
        friend arrayd operator*(double d, arrayd const &b){
            return b*d;
        }
        bool operator==(arrayd const &b) const{
            if(N==2){
                return array[0] == b[0] && array[1] == b[1];
            }
            if(N==3){
                return array[0] == b[0] && array[1] == b[1] && array[2] == b[2];
            }
            throw std::exception("Not implemented");
        }

        bool operator!=(arrayd const &b) const{
            return !(*this == b);
        }

        arrayd normalize() const{
            arrayd res(array);
            if(N==2){
                double dist = this->length();
                res[0] /= dist;
                res[1] /= dist;
                return res;
            }
            if(N==3){
                double dist = this->length();
                res[0] /= dist;
                res[1] /= dist;
                res[2] /= dist;
                return res;
            }
            throw std::exception("Not implemented");
        }
        
        double length() const{
            return sqrt(length2());
        }

        double length2() const{
            if(N==2){
                return array[0]*array[0] + array[1]*array[1];
            }
            if(N==3){
                return array[0]*array[0] + array[1]*array[1] + array[2]*array[2];
            }
            throw std::exception("Not implemented");
        }
        
        arrayd normal() const{
            arrayd res(array);
            double tmp = res[0]; res[0] = -res[1]; res[1] = tmp;
            return res;
        }

        double dot(arrayd const &a) const{
            if(N==2){
                return array[0]*a[0] + array[1]*a[1];
            }
            if(N==3){
                return array[0]*a[0] + array[1]*a[1] + array[2]*a[2];
            } 
            throw std::exception("Not implemented");
        }
        
        arrayd<3> cross(arrayd const &a) const{
            arrayd<3> res;
            if(N==2){
                res.z() = array[0]*a[1] - array[1]*a[0];
                return res;
            }
            throw std::exception("Not implemented");
        }
    };

    template<size_t N>
    struct SolverData{
        void reset(bool stream, bool grad, int size){
            if(stream){
                streamAcc.clear()
                streamError.clear();
            }
            if(grad){
                gradAcc.clear();
                gradError.clear();
            }
        }
        std::unordered_map<int,double> streamAcc;
        std::unordered_map<int,arrayd<N>> gradAcc;
        std::vector<double> streamError;
        std::vector<double> gradError;
    };

    template<size_t N>
    struct PDEData{
        void reset(int size){
            numSamples = 0;
            xAxis.clear();
        }
        std::mutex streamLock;
        int numSamples=0;
        std::unordered_map<int, arrayd<N>> computePoints;
        std::unordered_map<int, double> streamRef;
        std::unordered_map<int, arrayd<N>> gradRef;
        std::vector<int> xAxis;
    };
}

#endif