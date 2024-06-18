import matplotlib.pyplot as plt
import numpy as np

from smt.applications import EGO
from smt.applications.ego import Evaluator
from smt.surrogate_models import KRG, DesignSpace

def function_test_1d(x):
    # function xsinx
    import numpy as np

    x = np.reshape(x, (-1,))
    y = np.zeros(x.shape)
    y = (x - 3.5) * np.sin((x - 3.5) / (np.pi))
    return y.reshape((-1, 1))

n_iter = 3
n_parallel = 3
n_start = 50
xlimits = np.array([[0.0, 25.0]])

random_state = 42
design_space = DesignSpace(xlimits, random_state=random_state)
xdoe = np.atleast_2d([0, 7, 25]).T
n_doe = xdoe.size

class ParallelEvaluator(Evaluator):
    """
    Implement Evaluator interface using multiprocessing ThreadPool object (Python 3 only).
    """

    def run(self, fun, x):
        n_thread = 5
        # Caveat: import are made here due to SMT documentation building process
        from multiprocessing.pool import ThreadPool
        from sys import version_info

        import numpy as np

        if version_info.major == 2:
            return fun(x)
        # Python 3 only
        with ThreadPool(n_thread) as p:
            return np.array(
                [
                    y[0]
                    for y in p.map(
                        fun, [np.atleast_2d(x[i]) for i in range(len(x))]
                    )
                ]
            )

criterion = "EI"  #'EI' or 'SBO' or 'LCB'
qEI = "KBUB"  # "KB", "KBLB", "KBUB", "KBRand"
ego = EGO(
    n_iter=n_iter,
    criterion=criterion,
    xdoe=xdoe,
    surrogate=KRG(design_space=design_space, print_global=False),
    n_parallel=n_parallel,
    qEI=qEI,
    n_start=n_start,
    evaluator=ParallelEvaluator(),
    random_state=random_state,
)

x_opt, y_opt, _, x_data, y_data = ego.optimize(fun=function_test_1d)
print("Minimum in x={:.1f} with f(x)={:.1f}".format(x_opt.item(), y_opt.item()))

x_plot = np.atleast_2d(np.linspace(0, 25, 100)).T
y_plot = function_test_1d(x_plot)

fig = plt.figure(figsize=[10, 10])
for i in range(n_iter):
    k = n_doe + (i) * (n_parallel)
    x_data_k = x_data[0:k]
    y_data_k = y_data[0:k]
    x_data_sub = x_data_k.copy()
    y_data_sub = y_data_k.copy()
    for p in range(n_parallel):
        ego.gpr.set_training_values(x_data_sub, y_data_sub)
        ego.gpr.train()

        y_ei_plot = -ego.EI(x_plot)
        y_gp_plot = ego.gpr.predict_values(x_plot)
        y_gp_plot_var = ego.gpr.predict_variances(x_plot)

        x_data_sub = np.append(x_data_sub, x_data[k + p])
        y_KB = ego._get_virtual_point(np.atleast_2d(x_data[k + p]), y_data_sub)

        y_data_sub = np.append(y_data_sub, y_KB)

        ax = fig.add_subplot(n_iter, n_parallel, i * (n_parallel) + p + 1)
        ax1 = ax.twinx()
        (ei,) = ax1.plot(x_plot, y_ei_plot, color="red")

        (true_fun,) = ax.plot(x_plot, y_plot)
        (data,) = ax.plot(
            x_data_sub[: -1 - p],
            y_data_sub[: -1 - p],
            linestyle="",
            marker="o",
            color="orange",
        )
        (virt_data,) = ax.plot(
            x_data_sub[-p - 1 : -1],
            y_data_sub[-p - 1 : -1],
            linestyle="",
            marker="o",
            color="g",
        )

        (opt,) = ax.plot(
            x_data_sub[-1], y_data_sub[-1], linestyle="", marker="*", color="r"
        )
        (gp,) = ax.plot(x_plot, y_gp_plot, linestyle="--", color="g")
        sig_plus = y_gp_plot + 3.0 * np.sqrt(y_gp_plot_var)
        sig_moins = y_gp_plot - 3.0 * np.sqrt(y_gp_plot_var)
        un_gp = ax.fill_between(
            x_plot.T[0], sig_plus.T[0], sig_moins.T[0], alpha=0.3, color="g"
        )
        lines = [true_fun, data, gp, un_gp, opt, ei, virt_data]
        fig.suptitle("EGOp optimization of $f(x) = x \\sin{x}$")
        fig.subplots_adjust(hspace=0.4, wspace=0.4, top=0.8)
        ax.set_title("iteration {}.{}".format(i, p))
        fig.legend(
            lines,
            [
                "f(x)=xsin(x)",
                "Given data points",
                "Kriging prediction",
                "Kriging 99% confidence interval",
                "Next point to evaluate",
                "Expected improvment function",
                "Virtula data points",
            ],
        )
plt.show()