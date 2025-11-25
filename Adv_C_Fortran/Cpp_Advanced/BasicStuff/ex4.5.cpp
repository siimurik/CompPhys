#include <memory> // Requires C++11 or above
#include <iostream>

int main()
{
    std::shared_ptr<int> p_x(new int);
    std::cout<<"p_x use count: "<<p_x.use_count()<<"\n";
    *p_x = 5; // ’de-reference’ to alter contents
    
    // Create a weak_ptr pointing to the same resource
    std::weak_ptr<int> weak_p = p_x;
    std::cout << "weak_p use count: " << weak_p.use_count() << "\n";

    // Experiment 1: Access via weak_ptr BEFORE reset
    std::cout << "\n --- Before reset: ---\n";
    if (!weak_p.expired())
    {
        std::shared_ptr<int> temp = weak_p.lock();
        std::cout << "Value via weak_ptr (before reset): " << *temp << "\n";
        std::cout << "p_x use count with temp shared_ptr: " << p_x.use_count() << "\n";
    } else {
        std::cout << "weak_p is expired before reset.\n";
    }

    // Use this pointer elsewhere
    std::shared_ptr<int> p_y = p_x;
    std::cout<<"p_x use count: "<<p_x.use_count()<<"\n";
    p_y.reset();
    std::cout<<"p_x use count: "<<p_x.use_count()<<"\n";
    p_x.reset();
    std::cout<<"p_x use count: "<<p_x.use_count()<<"\n";

    // Experiment 2: Access via weak_ptr AFTER p_x.reset()
    p_x.reset(); // Now the resource is released
    std::cout << "\n --- After p_x.reset(): ---\n";
    std::cout << "p_x use count after reset: " << p_x.use_count() << "\n";

    if (!weak_p.expired())
    {
        std::shared_ptr<int> temp = weak_p.lock();
        std::cout << "Value via weak_ptr (after reset): " << *temp << "\n";
    } else 
    {
        std::cout << "Resource expired (after reset) - cannot acces value\n";
    }

    // Check what weak_p.lock() returns after reset
    std::shared_ptr<int> temp_after_reset = weak_p.lock();
    if (temp_after_reset)
    {
        std::cout << "Value via weak_ptr after reset: " << *temp_after_reset << "\n";
    }
    else
    {
        std::cout << "weak_p.lock() returned nullptr after reset.\n";
    }

    return 0;
}