#include <iostream>
using namespace std;

// Base class: Student
class Student {
private:
    double libraryFines; // Private to enforec non-negative constraint

public:
    string name;
    double tuitionFees;

    // Constructors
    Student() : name(""), libraryFines(0.0), tuitionFees(0.0) {}

    Student(string studentName)
        : name(studentName), libraryFines(0.0), tuitionFees(0.0) {}

    Student(string studentName, double fines, double tuition)
        : name(studentName), libraryFines(0.0), tuitionFees(tuition) {
            setLibraryFines(fines);
        }    

    // Setter for library fines (enforce non-negative values)
    void setLibraryFines(double fines) {
        if (fines >= 0.0) {
            libraryFines = fines;
        } else {
            cout << "Warning: Library fines must be non-negative. Setting to 0\n";
            libraryFines = 0.0;
        }
    }

    // Getter for library fines
    double getLibraryFines() const {
        return libraryFines;
    }

    // Virtual method to calculate total money owed
    virtual double getTotalMoneyOwed() const {
        return libraryFines + tuitionFees;
    }

    // Virtual destructor for proper cleanup in derived classes
    virtual ~Student() {}
};

// Derived class: GraduateStudent
class GraduateStudent : public Student {
private:
    bool isFullTime;  // true for full-time, false for part-time

public:
    // Constructors
    GraduateStudent() : Student(), isFullTime(true) {}

    GraduateStudent(string studentName, bool fullTime)
        : Student(studentName), isFullTime(fullTime) {}

    GraduateStudent(string studentName, double fines, bool fullTime)
        : Student(studentName, fines, 0.0), isFullTime(fullTime) {}

    // Setter and getter for enrollment status
    void setFullTimeStatus(bool fullTime) {
        isFullTime = fullTime;
    }

    bool getFullTimeStatus() const {
        return isFullTime;
    }

    // Override: Graduate students don't pay tuition fees
    double getTotalMoneyOwed() const override {
        return getLibraryFines(); // Only library fines, no tuition
    }
};

// Derived class: PhDStudent
class PhDStudent : public GraduateStudent {
public:
    // Constructors
    PhDStudent() : GraduateStudent() {}
    
    PhDStudent(string studentName, bool fullTime)
        : GraduateStudent(studentName, fullTime) {}
    
    PhDStudent(string studentName, double fines, bool fullTime)
        : GraduateStudent(studentName, fines, fullTime) {}

    // Override: Ph.D. students don't pay library fines or tuition
    double getTotalMoneyOwed() const override {
        return 0.0;  // Ph.D. students owe nothing
    }
};

int main() {
    
    cout << "=== University Student Management System ===\n\n";

    // Create an undergraduate student
    Student undergrad("Alice Johnson", 25.50, 15000.0);
    cout << "Undergraduate Student: " << undergrad.name << "\n";
    cout << "Library Fines: $" << undergrad.getLibraryFines() << "\n";
    cout << "Tuition Fees: $" << undergrad.tuitionFees << "\n";
    cout << "Total Owed: $" << undergrad.getTotalMoneyOwed() << "\n\n";

    // Create an undergraduate student
    GraduateStudent grad("Bob Smith", 15.00, true);
    cout << "Graduate Student: " << grad.name << "\n";
    cout << "Library Fines: $"   << grad.getLibraryFines() << "\n";
    cout << "Full-time: "        << (grad.getFullTimeStatus() ? "Yes" : "No") << "\n";
    cout << "Total Owed: $"      << grad.getTotalMoneyOwed() << "\n\n";

    // Create a Ph.D. student
    PhDStudent phd("Carol Davis", 50.00, true);
    cout << "Ph.D. Student: " << phd.name << endl;
    cout << "Library Fines (waived): $" << phd.getLibraryFines() << endl;
    cout << "Full-time: " << (phd.getFullTimeStatus() ? "Yes" : "No") << endl;
    cout << "Total Owed: $" << phd.getTotalMoneyOwed() << endl << endl;

    // Test non-negative enforcement
    cout << "=== Testing Non-negative Enforcement ===" << endl;
    Student testStudent("Test Student");
    testStudent.setLibraryFines(-10.0);  // Should show warning
    cout << "Library Fines after invalid set: $" << testStudent.getLibraryFines() << endl << endl;

    // Demonstrate polymorphism
    cout << "=== Demonstrating Polymorphism ===" << endl;
    Student* students[3];
    students[0] = &undergrad;
    students[1] = &grad;
    students[2] = &phd;

    for (int i = 0; i < 3; i++) {
        cout << "Student " << (i + 1) << " (" << students[i]->name 
             << ") owes: $" << students[i]->getTotalMoneyOwed() << endl;
    }

    return 0;
}