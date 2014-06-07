using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SwephNet
{
    /// <summary>
    /// Data
    /// </summary>
    public interface ISwephData
    {
        double EopTjdBeg { get; }
        double EopTjdBeg_horizons { get; }
        double EopTjdEnd { get; }
        double EopTjdEnd_add { get; }
        bool EopDpsiLoaded { get; }

        /// <summary>
        /// works for 100 years after 1962
        /// </summary>
        double[] Dpsi { get; }
        double[] Deps { get; }
    }
}
