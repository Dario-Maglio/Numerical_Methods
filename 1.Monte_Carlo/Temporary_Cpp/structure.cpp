#include <iostream>
#include <string>
using namespace std;

struct Person       //32 bytes   usually a CPU can acces to 64
{
    string name;    //24 bytes
    //else char[50] name;  and cin.get(variable, 50)
    int age;        //4 bytes
    float salary;   //4 bytes
};

int main()
{
    Person p1;

    cout << "Enter Full name: ";
    cin >> p1.name;
    cout << "Enter age: ";
    cin >> p1.age;
    cout << "Enter salary: ";
    cin >> p1.salary;

    cout << "\nDisplaying Information." << endl;
    cout << "Name: " << p1.name << endl;
    cout <<"Age: " << p1.age << endl;
    cout << "Salary: " << p1.salary << endl;

    return 0;
}
