from gaussian_dynamics import run_compact_v010_release_benchmark

out=run_compact_v010_release_benchmark()

print("v0.10 compact sensitivity/error budget")
print("----------------------------------------")
for key,value in out["error_budget"].items():
    print(f"{key}: {value}")

print("\nAcceptance result:")
print(out["acceptance"])

print(
    "\nThe sensitivity terms are not independent statistical errors and are "
    "not added in quadrature."
)
