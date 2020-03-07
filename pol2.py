# Polynomial Project.

# Start with the constant term followed by terms of higher order. Place zeros for any term that is not in use.


class Polynomial:

    def __init__(self, coeff):
        self.degree = len(coeff) - 1
        self.coeff = coeff

    def __gt__(self, other):
        if self.degree > other.degree:
            return True
        else:
            return False

    def __eq__(self, other):
        if self.degree == other.degree:
            return True
        else:
            return False

    def __lt__(self, other):
        if self.degree < other.degree:
            return True
        else:
            return False

    def __add__(self,other):

        if self.degree == other.degree:
            sum_coeff = [x + y for x, y in zip(self.coeff, other.coeff)]

        elif self.degree < other.degree:
            difference = other.degree - self.degree
            for i in range(difference):
                self.coeff.append(0)
            sum_coeff = [x + y for x, y in zip(self.coeff, other.coeff)]
        else:
            difference = self.degree - other.degree
            for i in range(difference):
                other.coeff.append(0)
            sum_coeff = [x + y for x, y in zip(self.coeff, other.coeff)]

        return Polynomial(sum_coeff)

    def __sub__(self,other):

        if self.degree == other.degree:
            sum_coeff = [x - y for x, y in zip(self.coeff, other.coeff)]

        elif self.degree < other.degree:
            difference = other.degree - self.degree
            for i in range(difference):
                self.coeff.append(0)
            sum_coeff = [x - y for x, y in zip(self.coeff, other.coeff)]
        else:
            difference = self.degree - other.degree
            for i in range(difference):
                other.coeff.append(0)
            sum_coeff = [x - y for x, y in zip(self.coeff, other.coeff)]

        return Polynomial(sum_coeff)

    def __repr__(self):
        return "Coefficients :{}".format(self.coeff)

    def differentiate(self):
        diff_coeff = []
        for i in range(self.degree + 1):
            diff_coeff.append(i)
        diff_pol = [a * b for a, b in zip(diff_coeff, self.coeff)]
        if len(diff_pol) == 0:
            print("Pol zero")
        else:
            diff_pol.pop(0)

        return Polynomial(diff_pol)

    def nr(self, iterations):

        y = self.differentiate()
        init_val = -99999999
        return_val = 0

        for iteration in range(1,iterations + 1):
            fx = 0
            fdx = 0
            for val in range(self.degree):
                fdx = fdx + (y.coeff[val] * (init_val ** val))

            for value in range(self.degree + 1):
                fx = fx + (self.coeff[value] * (init_val ** value))

            d = init_val - (fx / fdx)  # value of Newton-Raphson iteration
            init_val = d
            return_val = round(d, 5)

        init_val2 = 99999999
        return_val2 = 0

        for iteration in range(1, iterations + 1):
            fx = 0
            fdx = 0
            for val in range(self.degree):
                fdx = fdx + (y.coeff[val] * (init_val2 ** val))

            for value in range(self.degree + 1):
                fx = fx + (self.coeff[value] * (init_val2 ** value))

            d = init_val2 - (fx / fdx)  # value of Newton-Raphson iteration
            init_val2 = d
            return_val2 = round(d, 5)

        if return_val == return_val2:
            return 'tp is:', return_val
        else:
            return "tps are:", return_val, return_val2

    # Build tpsolve without for loops to increase efficiency. use map reduce filter
    # 1 list for coeff, 1 for




    def tpsolve(self):
        x = self.differentiate()

        if len(x.coeff) == 1:
            return "This function has no turning points"
        elif len(x.coeff) == 2:
            root = (-x.coeff[0])/x.coeff[1]
            return "Turning Point at:",root
        else:
            j = 0
            sol = []
            lf = 0.0001  # lf optimization current: 0.0001 # j,3 10M
            for val in range(-1000000, 1000000):  #(-10000000, 10000000)
                for i in range(self.degree):
                    j = j + (x.coeff[i]*((lf*val) ** i))
                if round(j,3) == 0:
                    sol.append(val*lf)
                j = 0
    # create a function that finds sign change over large intervals but does not miss any roots. NR??
            cleanup = []

            for item in range(1,len(sol)):
                if round((sol[item] - sol[item - 1]),2) == 0:
                    cleanup.append(item)

            cleanup.reverse()
            for item2 in cleanup:
                sol.pop(item2)

            return x.coeff, "tp(s) is(are):",sol


pol1 = Polynomial([0,0,0,3])
pol2 = Polynomial([2,4,6,8,10])
pol3 = Polynomial([1,2,3,4,5])
pol4 = Polynomial([1,-3,0,1])
pol5 = Polynomial([-1,1,-1,0,0,0,0,0,0,1,0,1])
pol6 = Polynomial([1,4,-5,9,6,23,14,5,6,1,1,10,1])
print(pol1 + pol2)

#print(pol2.differentiate())

print("with NR method")
print(pol1.differentiate().nr(100))
print(pol2.differentiate().nr(100))
print(pol3.differentiate().nr(100))
print(pol4.differentiate().nr(100))
print(pol5.differentiate().nr(100000))
print(pol6.differentiate().nr(10000))
print("with change of sign method")
print(pol1.tpsolve())
print(pol2.tpsolve())
print(pol3.tpsolve())
print(pol4.tpsolve())
print(pol5.tpsolve())

# RRE?

