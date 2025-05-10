import sys
from CPred_main import main

# Run the main function with adjusted arguments
if __name__ == "__main__":
    # The original code expects the first argument in the list
    sys.argv[0] = "CPred_main.py"  # This is just for consistency
    main()
