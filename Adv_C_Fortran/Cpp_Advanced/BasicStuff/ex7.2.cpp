#include <iostream>

// ============================================================================
// PART 1: ORIGINAL CODE
// ============================================================================
/*
class AbstractPerson
{
public:
    virtual void Print(){std::cerr<<"Never instantiate\n";}
};

class Mother : public AbstractPerson
{
public:
    virtual void Print(){std::cout<<"Mother\n";}
};

class Daughter : public Mother
{
public:
    void Print(){std::cout<<"Daughter\n";}
};

int main(int argc, char* argv[])
{
    std::cout << "=== PART 1: Original Code ===" << std::endl;
    AbstractPerson* p_mother = new Mother;
    AbstractPerson* p_daughter = new Daughter;
    p_mother->Print();      // Output: Mother
    p_daughter->Print();    // Output: Daughter
    delete p_mother;
    delete p_daughter;
    
    return 0;
}
*/

// ============================================================================
// PART 2: INVESTIGATING REMOVING 'public' KEYWORD FROM INHERITANCE
// ============================================================================
/*
// When you remove 'public' from inheritance, it becomes PRIVATE inheritance by default
// This means the base class members become private in the derived class

class AbstractPerson
{
public:
    virtual void Print(){std::cerr<<"Never instantiate\n";}
};

// Remove 'public' - inheritance becomes private by default
class Mother : AbstractPerson  // Private inheritance!
{
public:
    virtual void Print(){std::cout<<"Mother\n";}
};

class Daughter : public Mother
{
public:
    void Print(){std::cout<<"Daughter\n";}
};

int main(int argc, char* argv[])
{
    std::cout << "=== PART 2: Removing 'public' from Mother's inheritance ===" << std::endl;
    
    // THIS WILL NOT COMPILE!
    // Error: 'AbstractPerson' is an inaccessible base of 'Mother'
    // AbstractPerson* p_mother = new Mother;
    
    // The compiler cannot convert Mother* to AbstractPerson* because
    // the inheritance is private, making the base class inaccessible from outside.
    
    std::cout << "Code does not compile - private inheritance prevents base class pointer conversion" << std::endl;
    
    return 0;
}
*/

// ============================================================================
// PART 3: INVESTIGATING REMOVING 'virtual' KEYWORDS
// ============================================================================

// SCENARIO A: Remove 'virtual' from AbstractPerson::Print() (line 6)
/*
class AbstractPerson
{
public:
    void Print(){std::cerr<<"Never instantiate\n";}  // NOT virtual anymore
};

class Mother : public AbstractPerson
{
public:
    virtual void Print(){std::cout<<"Mother\n";}
};

class Daughter : public Mother
{
public:
    void Print(){std::cout<<"Daughter\n";}
};

int main(int argc, char* argv[])
{
    std::cout << "=== PART 3A: Remove 'virtual' from AbstractPerson::Print() ===" << std::endl;
    AbstractPerson* p_mother = new Mother;
    AbstractPerson* p_daughter = new Daughter;
    p_mother->Print();      // Output: "Never instantiate" (WRONG!)
    p_daughter->Print();    // Output: "Never instantiate" (WRONG!)
    delete p_mother;
    delete p_daughter;
    
    std::cout << "Without 'virtual' in base class, polymorphism doesn't work!" << std::endl;
    std::cout << "Static binding is used instead of dynamic binding." << std::endl;
    
    return 0;
}
*/

// SCENARIO B: Remove 'virtual' from Mother::Print() (line 12)
/*
class AbstractPerson
{
public:
    virtual void Print(){std::cerr<<"Never instantiate\n";}
};

class Mother : public AbstractPerson
{
public:
    void Print(){std::cout<<"Mother\n";}  // NOT virtual anymore
};

class Daughter : public Mother
{
public:
    void Print(){std::cout<<"Daughter\n";}
};

int main(int argc, char* argv[])
{
    std::cout << "=== PART 3B: Remove 'virtual' from Mother::Print() ===" << std::endl;
    AbstractPerson* p_mother = new Mother;
    AbstractPerson* p_daughter = new Daughter;
    p_mother->Print();      // Output: Mother (still works due to base class virtual)
    p_daughter->Print();    // Output: Daughter (still works)
    
    // The 'virtual' keyword in the base class is sufficient for polymorphism
    // Derived classes override virtual functions even without explicitly marking them virtual
    
    delete p_mother;
    delete p_daughter;
    
    return 0;
}
*/

// SCENARIO C: Add 'virtual' to Daughter::Print() (line 18)
/*
class AbstractPerson
{
public:
    virtual void Print(){std::cerr<<"Never instantiate\n";}
};

class Mother : public AbstractPerson
{
public:
    virtual void Print(){std::cout<<"Mother\n";}
};

class Daughter : public Mother
{
public:
    virtual void Print(){std::cout<<"Daughter\n";}  // Add 'virtual'
};

int main(int argc, char* argv[])
{
    std::cout << "=== PART 3C: Add 'virtual' to Daughter::Print() ===" << std::endl;
    AbstractPerson* p_mother = new Mother;
    AbstractPerson* p_daughter = new Daughter;
    p_mother->Print();      // Output: Mother
    p_daughter->Print();    // Output: Daughter
    delete p_mother;
    delete p_daughter;
    
    std::cout << "Adding 'virtual' is redundant but good practice for clarity." << std::endl;
    std::cout << "Modern C++ recommends using 'override' keyword instead." << std::endl;
    
    return 0;
}
*/

// ============================================================================
// PART 4: INSTANTIATING THE ABSTRACT CLASS
// ============================================================================
/*
class AbstractPerson
{
public:
    virtual void Print(){std::cerr<<"Never instantiate\n";}
};

class Mother : public AbstractPerson
{
public:
    virtual void Print(){std::cout<<"Mother\n";}
};

class Daughter : public Mother
{
public:
    void Print(){std::cout<<"Daughter\n";}
};

int main(int argc, char* argv[])
{
    std::cout << "=== PART 4: Instantiating AbstractPerson ===" << std::endl;
    
    // This WILL COMPILE and RUN (but shouldn't conceptually)
    AbstractPerson* p_abstract = new AbstractPerson;
    p_abstract->Print();    // Output: "Never instantiate"
    delete p_abstract;
    
    std::cout << "Problem: AbstractPerson can be instantiated!" << std::endl;
    std::cout << "This defeats the purpose of having an 'abstract' class." << std::endl;
    
    return 0;
}
*/

// ============================================================================
// PART 5: MAKING Print() PURE VIRTUAL
// ============================================================================
/*
class AbstractPerson
{
public:
    virtual void Print() = 0;  // Pure virtual function
};

class Mother : public AbstractPerson
{
public:
    virtual void Print(){std::cout<<"Mother\n";}
};

class Daughter : public Mother
{
public:
    void Print(){std::cout<<"Daughter\n";}
};

int main(int argc, char* argv[])
{
    std::cout << "=== PART 5: Pure Virtual Function ===" << std::endl;
    AbstractPerson* p_mother = new Mother;
    AbstractPerson* p_daughter = new Daughter;
    p_mother->Print();      // Output: Mother
    p_daughter->Print();    // Output: Daughter
    delete p_mother;
    delete p_daughter;
    
    std::cout << "Pure virtual function prevents instantiation of AbstractPerson." << std::endl;
    
    return 0;
}
*/

// ============================================================================
// PART 6: PURE VIRTUAL WITH REMOVED 'virtual' KEYWORDS
// ============================================================================
/*
// This scenario doesn't make sense: you cannot have a pure virtual function
// without the 'virtual' keyword. The syntax "= 0" requires 'virtual'.

// If you try:
class AbstractPerson
{
public:
    void Print() = 0;  // COMPILATION ERROR!
};

// Error: only virtual member functions can be declared pure
*/

// ============================================================================
// PART 7: TRYING TO INSTANTIATE PURE VIRTUAL CLASS
// ============================================================================

class AbstractPerson
{
public:
    virtual void Print() = 0;  // Pure virtual function
    virtual ~AbstractPerson() {}  // Virtual destructor for proper cleanup
};

class Mother : public AbstractPerson
{
public:
    virtual void Print(){std::cout<<"Mother\n";}
};

class Daughter : public Mother
{
public:
    void Print(){std::cout<<"Daughter\n";}
};

int main(int argc, char* argv[])
{
    std::cout << "=== PART 7: Attempting to Instantiate Pure Virtual Class ===" << std::endl;
    
    // THIS WILL NOT COMPILE!
    // Error: cannot instantiate abstract class
    // AbstractPerson* p_abstract = new AbstractPerson;
    
    std::cout << "Cannot instantiate AbstractPerson - it has pure virtual function(s)." << std::endl;
    std::cout << "This is the PROPER way to create an abstract class in C++." << std::endl << std::endl;
    
    // This still works fine:
    std::cout << "But derived classes with implementations work:" << std::endl;
    AbstractPerson* p_mother = new Mother;
    AbstractPerson* p_daughter = new Daughter;
    p_mother->Print();
    p_daughter->Print();
    delete p_mother;
    delete p_daughter;
    
    return 0;
}

/*
Key Findings:
Part 2 - Removing public keyword:

Result: Code won't compile
Reason: Without public, inheritance becomes private by default, making the base class inaccessible from outside, preventing pointer conversion

Part 3 - Removing virtual keywords:
Scenario A - Remove virtual from AbstractPerson::Print():

Result: Output is "Never instantiate" for both calls (WRONG!)
Reason: Without virtual in the base class, static binding is used instead of dynamic binding - polymorphism breaks completely

Scenario B - Remove virtual from Mother::Print():

Result: Still works correctly! Outputs "Mother" and "Daughter"
Reason: The virtual keyword in the base class is sufficient; derived classes automatically override virtual functions

Scenario C - Add virtual to Daughter::Print():

Result: No change in behavior
Reason: It's redundant but good practice for clarity (modern C++ uses override keyword instead)

Part 4 - Instantiating AbstractPerson:

Result: Compiles and runs, outputs "Never instantiate"
Problem: This defeats the purpose! The class isn't truly abstract

Part 5 - Pure Virtual Function:

Makes the class truly abstract
Cannot be instantiated
Derived classes must provide implementations

Part 6 - Pure Virtual without virtual:

Result: Compilation error
Reason: The syntax = 0 (pure virtual) requires the virtual keyword

Part 7 - Instantiating Pure Virtual Class:

Result: Compilation error - "cannot instantiate abstract class"
Conclusion: This is the PROPER way to create abstract classes in C++

Best Practice: Always use pure virtual functions (= 0) for abstract classes to prevent accidental instantiation and enforce the interface contract.
*/