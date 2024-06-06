# main.cs


using System;

public class ComplexNumber
{
    public double Real { get; }
    public double Imaginary { get; }

    public ComplexNumber(double real, double imaginary)
    {
        Real = real;
        Imaginary = imaginary;
    }

    public override string ToString()
    {
        return $"{Real} + {Imaginary}i";
    }
}

public class CubicEquationSolver
{
    public ComplexNumber[] Solve(double a, double b, double c)
    {
        double discriminant = Math.Pow(b, 2) - 4 * a * c;
        ComplexNumber[] roots;

        if (discriminant > 0)
        {
            roots = new ComplexNumber[3];
            double root1 = (-b + Math.Sqrt(discriminant)) / (2 * a);
            double root2 = (-b - Math.Sqrt(discriminant)) / (2 * a);
            roots[0] = new ComplexNumber(root1, 0);
            roots[1] = new ComplexNumber(root2, 0);
            roots[2] = new ComplexNumber(-(b / (2 * a)), 0);
        }
        else if (discriminant == 0)
        {
            roots = new ComplexNumber[2];
            double root = -b / (2 * a);
            roots[0] = new ComplexNumber(root, 0);
            roots[1] = new ComplexNumber(-(b / (2 * a)), 0);
        }
        else
        {
            roots = new ComplexNumber[3];
            double realPart = -b / (2 * a);
            double imaginaryPart = Math.Sqrt(-discriminant) / (2 * a);
            roots[0] = new ComplexNumber(realPart, imaginaryPart);
            roots[1] = new ComplexNumber(realPart, -imaginaryPart);
            roots[2] = new ComplexNumber(-(b / (2 * a)), 0);
        }

        return roots;
    }
}

public class Program
{
    public static void Main(string[] args)
    {
        Console.WriteLine("Введите коэффициенты кубического уравнения:");

        double[] coefficients;
        while (true)
        {
            coefficients = ReadCoefficients();
            if (coefficients.Length == 3)
                break;
            Console.WriteLine("Уравнение должно быть кубическим (3 коэффициента)");
        }

        CubicEquationSolver solver = new CubicEquationSolver();
        ComplexNumber[] roots = solver.Solve(coefficients[0], coefficients[1], coefficients[2]);

        Console.WriteLine("Корни кубического уравнения:");
        foreach (var root in roots)
        {
            Console.WriteLine(root);
        }
    }

    private static double[] ReadCoefficients()
    {
        double[] coefficients = new double[3];
        string[] coefficientNames = { "a", "b", "c" };

        for (int i = 0; i < coefficients.Length; i++)
        {
            while (true)
            {
                Console.Write($"Введите коэффициент {coefficientNames[i]}: ");
                if (double.TryParse(Console.ReadLine(), out double coefficient))
                {
                    coefficients[i] = coefficient;
                    break;
                }
                Console.WriteLine("Некорректный ввод. Пожалуйста, введите число.");
            }
        }

        return coefficients;
    }
}

# equation.cs
using System;
using System.Linq;
using System.Numerics;

public interface IPolynomial
{
    int Dimension { get; }
    double[] Coefficients { get; }
    Complex[] FindRoots();
}

public interface ISolutionStrategy
{
    Complex[] Solve(double[] coefficients);
}

public class Polynomial : IPolynomial
{
    private readonly double[] _coefficients;
    private readonly ISolutionStrategy _solutionStrategy;

    public Polynomial(double[] coefficients, ISolutionStrategy solutionStrategy)
    {
        _coefficients = coefficients;
        _solutionStrategy = solutionStrategy;
    }

    public int Dimension => _coefficients.Length;
    public double[] Coefficients => (double[])_coefficients.Clone();

    public Complex[] FindRoots()
    {
        double[] coefficients = _coefficients.Select(x => (double)x).ToArray();
        return _solutionStrategy.Solve(coefficients);
    }
}

public class LinearStrategy : ISolutionStrategy
{
    public Complex[] Solve(double[] coefficients)
    {
        if (coefficients.Length != 2)
        {
            throw new ArgumentException("Неверное количество коэффициентов для линейного уравнения.");
        }
        double a = coefficients[0];
        double b = coefficients[1];
        if (a == 0)
        {
            throw new InvalidOperationException("Комплексные корни отсутствуют");
        }
        return new[] { new Complex(-b / a, 0) };
    }
}

public class QuadraticStrategy : ISolutionStrategy
{
    public Complex[] Solve(double[] coefficients)
    {
        if (coefficients.Length != 3)
        {
            throw new ArgumentException("Неверное количество коэффициентов для квадратного уравнения.");
        }
        double a = coefficients[0], b = coefficients[1], c = coefficients[2];
        double discriminant = b * b - 4 * a * c;
        if (discriminant >= 0)
        {
            double x1 = (-b + Math.Sqrt(discriminant)) / (2 * a);
            double x2 = (-b - Math.Sqrt(discriminant)) / (2 * a);
            return new[] { new Complex(x1, 0), new Complex(x2, 0) };
        }
        else
        {
            double realPart = -b / (2 * a);
            double imaginaryPart = Math.Sqrt(-discriminant) / (2 * a);
            return new[] { new Complex(realPart, imaginaryPart), new Complex(realPart, -imaginaryPart) };
        }
    }
}

public static class Equations
{
    public static Polynomial CreateEquation(double[] coefficients)
    {
        var trimmedCoefficients = TrimCoefficients(coefficients);
        var strategy = ChooseStrategy(trimmedCoefficients);
        return new Polynomial(trimmedCoefficients, strategy);
    }

    private static double[] TrimCoefficients(double[] coefficients)
    {
        int nonZeroIndex = Array.FindLastIndex(coefficients, x => x != 0);
        if (nonZeroIndex == -1) return new double[] { 0 };
        return coefficients.Take(nonZeroIndex + 1).ToArray();
    }

    private static ISolutionStrategy ChooseStrategy(double[] coefficients)
    {
        switch (coefficients.Length - 1)
        {
            case 0:
                throw new InvalidOperationException("Корней бесконечно много");
            case 1:
                return new LinearStrategy();
            case 2:
                return new QuadraticStrategy();
            default:
                throw new InvalidOperationException("Неизвестный тип уравнения или слишком много коэффициентов");
        }
    }
}

public class Program
{
    public static void Main(string[] args)
    {

        double[] coefficients = { 1, -3, 2 }; 
        var equation = Equations.CreateEquation(coefficients);
        var roots = equation.FindRoots();
        Console.WriteLine("Roots:");
        foreach (var root in roots)
        {
            Console.WriteLine(root);
        }
    }
}

# complex.cs
using System;

public class Complex
{

    public double Re; 
    public double Im; 

    public Complex(double re, double im)
    {
        Re = re;
        Im = im;
    }

    public Complex Add(Complex other)
    {
        return new Complex(this.Re + other.Re, this.Im + other.Im);
    }

    public Complex Subtract(Complex other)
    {
        return new Complex(this.Re - other.Re, this.Im - other.Im);
    }

    public Complex Multiply(Complex other)
    {
        double newRe = this.Re * other.Re - this.Im * other.Im;
        double newIm = this.Re * other.Im + this.Im * other.Re;
        return new Complex(newRe, newIm);
    }

    public Complex Divide(Complex other)
    {
        if (other.Re == 0 && other.Im == 0)
        {
            throw new Exception("Деление на ноль.");
        }

        double denom = other.Re * other.Re + other.Im * other.Im;
        double newRe = (this.Re * other.Re + this.Im * other.Im) / denom;
        double newIm = (this.Im * other.Re - this.Re * other.Im) / denom;
        return new Complex(newRe, newIm);
    }

    public override string ToString()
    {
        return String.Format("{0} + {1}i", Re, Im);
    }
}

class Program
{
    static void Main()
    {
        Complex num1 = new Complex(2, 3);
        Complex num2 = new Complex(1, 1);

        Complex resultAdd = num1.Add(num2);
        Console.WriteLine("Сложение: " + resultAdd.ToString());

        Complex resultSubtract = num1.Subtract(num2);
        Console.WriteLine("Вычитание: " + resultSubtract.ToString());

        Complex resultMultiply = num1.Multiply(num2);
        Console.WriteLine("Умножение: " + resultMultiply.ToString());

        Complex resultDivide = num1.Divide(num2);
        Console.WriteLine("Деление: " + resultDivide.ToString());
    }
}
# strategy.cs
using System;

public interface IPolynomialSolver
{
    double[] Solve(double[] coefficients);
}

public static class Strategies
{
    public const string Linear = "Линейное уравнение";
    public const string Quadratic = "Квадратное уравнение";

    public static IPolynomialSolver GetSolver(string strategy)
    {
        switch (strategy)
        {
            case Linear:
                return new LinearEquationSolver();
            case Quadratic:
                return new QuadraticEquationSolver();
            default:
                throw new ArgumentException("Неизвестная стратегия решения уравнения");
        }
    }
}

public class LinearEquationSolver : IPolynomialSolver
{
    public double[] Solve(double[] coefficients)
    {
        double a = coefficients[0];
        double b = coefficients[1];

        if (a == 0)
        {
            throw new ArgumentException("Коэффициент 'a' не может быть равен нулю");
        }

        return new double[] { -b / a };
    }
}

public class QuadraticEquationSolver : IPolynomialSolver
{
    public double[] Solve(double[] coefficients)
    {
        double a = coefficients[0];
        double b = coefficients[1];
        double c = coefficients[2];

        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0)
        {
            throw new ArgumentException("Уравнение не имеет действительных корней");
        }

        double sqrtDiscriminant = Math.Sqrt(discriminant);
        double x1 = (-b + sqrtDiscriminant) / (2 * a);
        double x2 = (-b - sqrtDiscriminant) / (2 * a);

        return new double[] { x1, x2 };
    }
}

public class Program
{
    public static void Main(string[] args)
    {
        Console.WriteLine("Введите коэффициенты полинома:");

        double[] coefficients;
        while (true)
        {
            coefficients = ReadCoefficients();
            if (coefficients.Length >= 2)
                break;
            Console.WriteLine("Уравнение должно быть как минимум линейным (2 коэффициента)");
        }

        Console.WriteLine("Выберите стратегию решения:");
        Console.WriteLine($"1. {Strategies.Linear}");
        Console.WriteLine($"2. {Strategies.Quadratic}");

        string choice = Console.ReadLine();
        IPolynomialSolver solver;

        try
        {
            solver = GetSolver(choice);
        }
        catch (ArgumentException ex)
        {
            Console.WriteLine(ex.Message);
            return;
        }

        try
        {
            double[] roots = solver.Solve(coefficients);
            Console.WriteLine("Корни уравнения:");
            foreach (var root in roots)
            {
                Console.WriteLine(root);
            }
        }
        catch (ArgumentException ex)
        {
            Console.WriteLine(ex.Message);
        }
    }

    private static double[] ReadCoefficients()
    {
        double[] coefficients = new double[3];
        string[] coefficientNames = { "a", "b", "c" };

        for (int i = 0; i < coefficients.Length; i++)
        {
            while (true)
            {
                Console.Write($"Введите коэффициент {coefficientNames[i]}: ");
                if (double.TryParse(Console.ReadLine(), out double coefficient))
                {
                    coefficients[i] = coefficient;
                    break;
                }
                Console.WriteLine("Некорректный ввод. Пожалуйста, введите число.");
            }
        }

        return coefficients;
    }

    private static IPolynomialSolver GetSolver(string choice)
    {
        switch (choice)
        {
            case "1":
                return Strategies.GetSolver(Strategies.Linear);
            case "2":
                return Strategies.GetSolver(Strategies.Quadratic);
            default:
                throw new ArgumentException("Неверный выбор стратегии");
        }
    }
}
